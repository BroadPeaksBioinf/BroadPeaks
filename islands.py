#!/usr/bin/env/ python

import numpy
import scipy
import pysam
import logging

"""
This code repeats, maybe function like this may help.

def yield_all_reads_in_chromosome(chromosomes_info):
     for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        current_chromosome_name = chromosome[0]
        # currentChromosomeSize = int(chromosome[1])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
        yield all_reads_in_chromosome
"""
# Makes a simple list of windows, where each window is a list [WindowStart, ReadsInWindow].
# Chromosomes are separated by [-1,-1] window
# Sequences of ineligible windows longer than GAP+1 are not stored


def make_windows_list(bam_path, chromosomes_info, l0, window_size, gap, unique_reads_count, normalization_coef):
    logging.info("Making eligible windows of {} bp with allowed gap_size {} bp".format(window_size, window_size*gap))
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    window_list = []
    # logging.info("chromosome_name, chromosome_size, total_number_of_eligible_windows_on_chromosome")
    for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        previous_read_strand = 0
        current_chromosome_name = chromosome[0]
        current_chromosome_size = int(chromosome[1])
        # print([current_chromosome_name, current_chromosome_size, len(window_list)])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
        gap_count = 0
        window_start = 0
        window_reads_count = 0
        # just for precise number of windows counting
        chr_window = 0
        i = 0
        for read in all_reads_in_chromosome:
            read_str = str(read)
            # read strand: 0 = +         16 = -
            read_strand = ([int(s) for s in read_str.split() if s.isdigit()][0])
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            # filtering redundant reads
            if (beginning_of_the_read != beginning_of_the_previous_read) or (read_strand != previous_read_strand):
                beginning_of_the_previous_read = beginning_of_the_read
                previous_read_strand = read_strand
                # Ground state: gap_count <= GAP
                gap_flag = True
                while True:
                    # if read in window
                    if window_start <= beginning_of_the_read and beginning_of_the_read < window_start + window_size:
                        window_reads_count += 1
                        break
                    # if read before window: NEVER ENTERING THIS CONDITION
                    elif beginning_of_the_read < window_start:
                        break
                    else:
                        window_reads_count = int(float(window_reads_count)*normalization_coef)
                        if window_reads_count < l0:
                            gap_count += 1
                        else:
                            gap_flag = False
                            gap_count = 0
                        # * unique_reads_count/1000000 is for normalization per million reads
                        # now we are able to compare control and sample

                        # NEED TO CHANGE LAMBDA AND N
                        window_list.append([window_start, window_reads_count])
                        chr_window += 1
                        # / (unique_reads_count/1000000)])
                        # If we have a g+1 sized GAP, go and delete last g windows
                        if gap_count > gap or gap_flag:
                            gap_flag = True
                            while gap_count > 0:
                                window_list.pop()
                                gap_count -= 1
                                chr_window -= 1
                        window_start += window_size
                        window_reads_count = 0
        # Next chromosome marker just in case
        window_list.append([-1, -1])
        i += 1
        logging.info("On {} there are {} eligible windows".format(current_chromosome_name, chr_window))
    # candidate windows == eligible + ineligible(below gap)
    logging.info("\nThere are {} candidate windows".format(len(window_list) - i))
    window_list.append([1, 1])
    bamfile.close()
    return window_list




def modify_window_list_based_on_control(control_path, chromosomes_info, l0, window_size, gap, unique_reads_count, control_unique_reads_count, window_list_wo_control):
    logging.info("Making window list based on control")
    bamfile = pysam.AlignmentFile(control_path, 'rb')
    NORMALIZATION_CONSTANT = float(unique_reads_count)/float(control_unique_reads_count)
    window_list = []
    i = 0
    # logging.info("chromosome_name, chromosome_size, total_number_of_eligible_windows_on_chromosome")
    for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        previous_read_strand = 0
        current_chromosome_name = chromosome[0]
        current_chromosome_size = int(chromosome[1])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
        window_start = window_list_wo_control[i][0]
        control_window_reads_count = 0

        for read in all_reads_in_chromosome:
            read_str = str(read)
            # read strand: 0 = +         16 = -
            read_strand = ([int(s) for s in read_str.split() if s.isdigit()][0])
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            # filtering redundant reads
            if (beginning_of_the_read != beginning_of_the_previous_read) or (read_strand != previous_read_strand):
                beginning_of_the_previous_read = beginning_of_the_read
                previous_read_strand = read_strand

                if window_start <= beginning_of_the_read and beginning_of_the_read < window_start + window_size:
                    control_window_reads_count += 1
                    break
                elif beginning_of_the_read < window_start:
                    break
                else:
                    window_reads_count = window_list_wo_control[i][1] - int(control_window_reads_count*NORMALIZATION_CONSTANT)
                    window_list.append([window_start, window_reads_count])
                    i += 1
                    window_start = window_list_wo_control[i][0]
                    #next chromosome check
                    if window_start == -1:
                        window_list.append([-1,-1])
                        window_start = window_list_wo_control[i][0]
                        i += 1
                        break
        #go and do next chrom

    #prune new window list of gaps
    gap_count = 0
    gap_flag = True
    window_list_new = []
    for window in window_list:
        window_list_new.append(window)
        if window[1] < l0:
            gap_count += 1
        else:
            gap_flag = False
        if gap_count > gap or gap_flag:
            gap_flag = True
            while gap_count > 0:
                window_list_new.pop()
                gap_count -= 1

    return window_list_new




def calculate_window_score(reads_in_window, lambdaa, l0):
    # sometimes 0 and therefore inf in -log  is generated
    if reads_in_window >= l0:
        temp = scipy.stats.poisson.pmf(reads_in_window, lambdaa)
        if temp < 1e-320:
            window_score = 1000
        else:
            window_score = -numpy.log(temp)
    else:
        window_score = 0
    return window_score




def make_islands_list(window_list, lambdaa, window_size, l0, chromosomes_info, island_score_threshold):
    chromosome_counter = 0
    current_chromosome_name = chromosomes_info[chromosome_counter][0]
    islands_list = []
    island_score = 0
    island_number_of_reads = 0
    island_number_of_gaps = 0
    window_start = window_list[0][0] - window_size
    island_start = window_list[0][0]
    zero_score_islands = 0
    for i, window in enumerate(window_list):
        # i == # in list, window == [window_start_position, number_of_reads_per_window]
        window_start_new = window[0]
        # what is window[1]?
        number_of_reads = window[1]

        # New chromosome check: [-1  -1] window separates chomosomes.
        if window_start_new == -1:
            window_start = window_list[i + 1][0] - window_size
            chromosome_counter += 1
            # switch the chromosome name to next one
            if chromosome_counter < len(chromosomes_info):
                current_chromosome_name = chromosomes_info[chromosome_counter][0]
        else:
            # if the next window belongs to the next island:
            if window_start_new != window_start + window_size:
                # Special case for the one-window island:
                if window_start == island_start:
                    island_score = calculate_window_score(number_of_reads, lambdaa, l0)
                    island_number_of_reads = number_of_reads
                island_length = (window_start + window_size - island_start)/window_size
                islands_list.append([current_chromosome_name, island_start,
                                    window_start + window_size, island_score, island_number_of_reads, island_length - island_number_of_gaps, island_number_of_gaps])
                island_score = 0
                island_number_of_reads = 0
                island_number_of_gaps = 0
                island_start = window_start_new
            else:
                # Set score
                window_score = calculate_window_score(number_of_reads, lambdaa, l0)
                if number_of_reads<l0:
                    island_number_of_gaps += 1
                island_score += window_score
                island_number_of_reads += number_of_reads
            window_start = window_start_new
    logging.info("There are {} islands found".format(len(islands_list)))
    print(zero_score_islands)
    return islands_list



#island_list = [[0,1,5], [0,10,12], [0,14,15],[0,16,17], [0,25,27], [0,100,101], [0,200,201],[0,1000,1001],[0,1000000,10000001]]
#island_list_2 = [[0,1,2], [0,6,7], [0,8,9], [0,11,12], [0,14,15], [0, 99,102],[0,500,501],[0,1002,1003],[0,1005,1006],[0,1007,1008]]

def find_unintersected_islands(island_list, island_list_2):

    final_islands = []
    i = 0
    second_island_beginning = island_list_2[i][1]
    second_island_end = island_list_2[i][2]
    # whole second island before first island flag
    flag = 0

    for (j,island) in enumerate(island_list):
        intersection_flag = 0
        first_island_beginning = island[1]
        first_island_end = island[2]
        #print(first_island_beginning,first_island_end,second_island_beginning,second_island_end)

        if i>=len(island_list_2)-1:
            final_islands.append(island)
            continue

        if second_island_beginning > first_island_end:
            final_islands.append(island)
            continue


        while True:
            if i>=len(island_list_2):
                break
            if second_island_beginning > first_island_end:
                flag = 1
                break
            second_island_beginning = island_list_2[i][1]
            second_island_end = island_list_2[i][2]
            #print(first_island_beginning,first_island_end,second_island_beginning,second_island_end)

            if (second_island_end < first_island_beginning):
                flag = 0
                i +=1
            # intersection condition
            elif((second_island_beginning>=first_island_beginning) and (second_island_beginning<=first_island_end)) or\
                    ((second_island_end>=first_island_beginning) and (second_island_end<=first_island_end)) or \
                    ((second_island_beginning<first_island_beginning) and (second_island_end>first_island_end)):
                intersection_flag = 1
                i+=1

        if flag == 1 and intersection_flag == 0:
            final_islands.append(island)
            flag = 0
        intersection_flag = 0

        if i>=len(island_list_2)-1:
            final_islands.append(island)
            continue

    return (final_islands)

#final_i = find_unintersected_islands(island_list, island_list_2)
#print(final_i)

