#!/usr/bin/env/ python

import numpy
import scipy
import pysam


# Makes a simple list of windows, where each window is a list [WindowStart, ReadsInWindow].
# Chromosomes are separated by [-1,-1] window
# Sequences of ineligible windows longer than GAP+1 are not stored
def make_windows_list(bam_path, chromosomes_info, l0, window_size, gap, unique_reads_count):
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    window_list = []
    for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        current_chromosome_name = chromosome[0]
        current_chromosome_size = int(chromosome[1])
        # print([current_chromosome_name, current_chromosome_size, len(window_list)])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)

        gap_count = 0
        window_start = 0
        window_reads_count = 0
        for read in all_reads_in_chromosome:
            read_str = str(read)
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            # filtering redundant reads
            if beginning_of_the_read != beginning_of_the_previous_read:
                beginning_of_the_previous_read = beginning_of_the_read
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
                        if window_reads_count < l0:
                            gap_count += 1
                        else:
                            gap_flag = False
                            gap_count = 0
                        # * unique_reads_count/1000000 is for normalization per million reads
                        # now we are able to compare control and sample

                        # NEED TO CHANGE LAMBDA AND N
                        window_list.append([window_start, window_reads_count])
                        # / (unique_reads_count/1000000)])
                        # If we have a g+1 sized GAP, go and delete last g windows
                        if gap_count > gap or gap_flag:
                            gap_flag = True
                            while gap_count > 0:
                                window_list.pop()
                                gap_count -= 1
                        window_start += window_size
                        window_reads_count = 0
        # Next chromosome marker just in case
        window_list.append([-1, -1])
        print([current_chromosome_name, current_chromosome_size, len(window_list)], "READY")
    window_list.append([1, 1])
    bamfile.close()
    return window_list






def make_islands_list(window_list, lambdaa, window_size, l0, chromosomes_info, island_score_threshold):
    chromosome_counter = 0
    current_chromosome_name = chromosomes_info[chromosome_counter][0]
    islands_list = []
    island_score = 0
    window_start = window_list[0][0] - window_size
    island_start = window_list[0][0]

    for i, window in enumerate(window_list):
        # i == # in list, window == [window_start_position, number_of_reads_per_window]
        window_start_new = window[0]
        # what is window[1]?
        number_of_reads = window[1]

        # New chromosome check
        if window_start_new == -1:
            print (current_chromosome_name + " done")
            # move to the next-previous? window
            window_start = window_list[i + 1][0] - window_size
            chromosome_counter += 1

            if chromosome_counter < len(chromosomes_info):
                current_chromosome_name = chromosomes_info[chromosome_counter][0]
                # print ("start " + current_chromosome_name)
        else:
            # + window_size * gap?
            # if window_start_new < window_start * (gap + 1) ?
            if window_start_new != window_start + window_size:
                # A bug here: loads of 0-score islands are generated
                if island_score >= island_score_threshold:
                    islands_list.append([current_chromosome_name, island_start,
                                         window_start + window_size, int(island_score)])
                island_score = 0
                island_start = window_start_new
            else:
                # Check eligibility
                if number_of_reads >= l0:
                    # sometimes 0 and therefore inf in -log  is generated
                    temp = scipy.stats.poisson.pmf(number_of_reads, lambdaa)
                    # what for temp == 0?
                    if temp == 0:
                        window_score = 10
                    else:
                        window_score = -numpy.log(temp)
                else:
                    window_score = 0
                island_score += window_score
            window_start = window_start_new

    return islands_list
