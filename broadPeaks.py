#!/usr/bin/env/ python

import sys
import numpy
import pysam
import logging
import os.path
import time
import scipy
import scipy.stats


def input_and_index(bam_path):
    # Check if there is an BAM file at entered dir
    if not os.path.isfile(bam_path):
        logging.error("BAM file does not exist")

    # Check if there is an index file, create one if there isn't
    if not os.path.isfile(bam_path + ".bai"):
        pysam.index(bam_path)
        logging.info('No index was found, new index was generated')

    # Take chromosome data from BAM index:
    # (ref.seq. name, ref.seq. length, number of mapped reads and number of unmapped reads)
    chromosomes_info = []
    for chr in pysam.idxstats(bam_path):
        chromosomes_info.append(chr.split("\t"))
    # Last line is unmapped reads, we don't need them
    chromosomes_info.pop()
    # print(chromosomes_info)

    bamfile = pysam.AlignmentFile(bam_path, "rb")

    return (bamfile, chromosomes_info)


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


def count_unique_reads(bamfile, chromosomes_info):
    total_unique_reads_count = 0
    for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        current_chromosome_name = chromosome[0]
        # currentChromosomeSize = int(chromosome[1])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
        for read in all_reads_in_chromosome:
            read_str = str(read)
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            if beginning_of_the_read != beginning_of_the_previous_read:
                beginning_of_the_previous_read = beginning_of_the_read
                total_unique_reads_count += 1
    # print("Unique reads counted")
    return total_unique_reads_count
    # normalising_coefficient = total_unique_reads_count / 1000000
    # it can help to calculate experiments with "control"
    # read_coverage has to be multiplied on normalising_coefficient


# Makes a simple list of windows, where each window is a list [WindowStart, ReadsInWindow].
# Chromosomes are separated by [-1,-1] window
# Sequences of ineligible windows longer than GAP+1 are not stored


def make_windows_list(bamfile, chromosomes_info, l0, window_size, gap):
    window_list = []
    for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        current_chromosome_name = chromosome[0]
        current_chromosome_size = int(chromosome[1])
        # print([current_chromosome_name, current_chromosome_size, len(window_list)])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)

        gap_count = 0
        i = 0
        window_reads_count = 0
        for read in all_reads_in_chromosome:
            read_str = str(read)
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            if beginning_of_the_read != beginning_of_the_previous_read:
                beginning_of_the_previous_read = beginning_of_the_read
                # Ground state: gap_count <= gap
                gap_flag = 1
                while True:
                    # condition seems to be equal, but not
                    # if i < i + window_size
                    if i <= beginning_of_the_read and beginning_of_the_read < i + window_size:
                        window_reads_count += 1
                        break
                    elif beginning_of_the_read < i:
                        break
                    else:
                        if window_reads_count < l0:
                            gap_count += 1
                        else:
                            gap_flag = 0
                        window_list.append([i, window_reads_count])
                        # If we have a g+1 sized GAP, go and delete last g windows
                        if gap_count > gap or gap_flag == 1:
                            gap_flag = 1
                            while gap_count > 0:
                                window_list.pop()
                                gap_count -= 1
                        i += window_size
                        window_reads_count = 0

        # Next chromosome marker just in case
        window_list.append([-1, -1])
        print([current_chromosome_name, current_chromosome_size, len(window_list)], "READY")
    window_list.append([1, 1])
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
                    if temp == 0:
                        window_score = 10
                    else:
                        window_score = -numpy.log(temp)
                else:
                    window_score = 0
                island_score = island_score + window_score
            window_start = window_start_new

    return islands_list


startTime = time.time()

# Input arguments: path to BAM file
# bamPath = sys.argv[1]
# WINDOW_SIZE = sys.argv[2]
# GAP = sys.argv[3]
# GAP = int(float(GAP)/float(WINDOW_SIZE))
# bamPath = "/home/dima/BAMfiles/Bernstein_H1_hESC_CTCF.bam"
# "/home/dima/BAMfiles/h3k4me3_rep1.bam"
bamPath = "/home/yegor/Alex_project/H3K4me3.bam"
WINDOW_SIZE = 200
p0 = 0.01
GAP = 1
ISLAND_SCORE_THRESHOLD = 100

# Log file
logging.basicConfig(filename=(os.path.dirname(bamPath) + '/SICER_log.log'), level=logging.DEBUG)

bamfile, chromosomes_info = input_and_index(bamPath)
print("Counting unique reads")
total_unique_reads_count = count_unique_reads(bamfile, chromosomes_info)

# Effective genome length
L = 0.77 * sum(int(row[1]) for row in chromosomes_info)
# Lambda for poisson dist
lambdaa = float(WINDOW_SIZE) * float(total_unique_reads_count) / float(L)
# Minimum #reads in a window for eligibility
# Formula (1), finding l0
l0 = scipy.stats.poisson.ppf(1 - p0, lambdaa)
# print(L, total_unique_reads_count, lambdaa, l0)


print("Finished counting reads, now making window list")
windowList = make_windows_list(bamfile, chromosomes_info, l0, WINDOW_SIZE, GAP)
print("Finished window list, now making island list")
island_list = make_islands_list(windowList, lambdaa, WINDOW_SIZE, l0, chromosomes_info,
                                ISLAND_SCORE_THRESHOLD)

# print(len(windowList), sys.getsizeof(windowList)/1024)
# print(len(island_list), sys.getsizeof(island_list)/1024)

f = open(bamPath[:-4] + '_peaks.bed', 'wb')
for island in island_list:
    islandString = str(island[0]) + "\t" + str(island[1]) + "\t" + str(island[2]) + "\n"
    f.write(islandString)
# print(island_list[11104])

bamfile.close()

print("Finished. Elapsed time, minutes: " + str((time.time() - startTime) / 60))
