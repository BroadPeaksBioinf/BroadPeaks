#!/usr/bin/env/ python

import sys
import numpy
import pysam
import logging
import os.path
import time
import scipy
import scipy.stats
import argparse


def input_and_index(bam_path):
    # Check if there is a file at entered dir
    if not os.path.isfile(bam_path):
        logging.error("No BAM file specified in input '{}' or there is no such a file".format(bam_path))
        sys.exit("`{}` is not a path to BAM file. \n More information in `{}`".format(bam_path, LOG_filename))

    if bam_path[-4:] != '.bam':
        logging.error("`{}` is other file type, not BAM. This tool works only with BAM-files as input".
                      format(bam_path))
        sys.exit("`{}` is not a BAM file. \n More information in `{}`".format(bam_path, LOG_filename))

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
                # Ground state: gap_count <= GAP
                gap_flag = 1
                while True:
                    # condition seems to be equal, but not
                    # if i < i + WINDOW_SIZE
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

parser = argparse.ArgumentParser(description="Tool for ChiP-seq analysis to find broad peaks")

parser.add_argument('infile', help="Path to `input.bam` file", type=str)
parser.add_argument('-w', dest='window_size', help="Window size (bp).  DEFAULT: 200", type=int, default=200)
parser.add_argument('-g', dest='gap', help="Gap size shows how many windows could be skipped. DEFAULT: 1",
                    type=int, default=1, choices=[1, 2, 3])
parser.add_argument('-p', dest='p_value', help="p-value; has to be in range(0.0, 1.0). DEFAULT: 0.01",
                    type=float, default=0.01)
parser.add_argument('-t', dest='threshold', help="Island score threshold. DEFAULT: 100", type=int, default=100)
parser.add_argument('-o', dest='outfile', help="Path to `output.bed` file. DEFAULT: path to `input.bam` file; "
                                               "output file name is `input_peaks.bed`", type=file)
parser.add_argument('-e', help="Proportion of effective genome length; has to be in range(0.0, 1.0) DEFAULT: 0.77",
                    type=float, default=0.77)
parser.add_argument('-c', dest='control', help="Path to `control.bam` file. DEFAULT: no control file",
                    type=str)
parser.add_argument('-ref', dest='reference genome', help="Reference genome.  DEFAULT: 'hg19'",
                    type=str, default='hg19')

args = parser.parse_args(["/home/yegor/Alex_project/H3K4me3.bam"])

# Input arguments: see broadPeaks.py -h or --help
# bamPath = "/home/dima/BAMfiles/Bernstein_H1_hESC_CTCF.bam"
# "/home/dima/BAMfiles/h3k4me3_rep1.bam"
bamPath = args.infile
WINDOW_SIZE = args.window_size
GAP = args.gap
p0 = args.p_value
if p0 <= 0 or p0 >= 1:
    logging.error("p-value has to be in range(0.0, 1.0)")
if args.e <= 0 or args.e > 1:
    logging.error("proportion of effective genome length has to be in range(0.0, 1.0)")
ISLAND_SCORE_THRESHOLD = args.threshold
outfile = args.outfile
if not outfile:
    outfile = bamPath[:-4] + '_peaks.bed'

controlPath = args.control
LOG_filename = 'SICER_log.log'

# Log file
logging.basicConfig(filename=(os.path.dirname(bamPath) + '/' + LOG_filename), level=logging.DEBUG,
                    format='%(asctime)s : %(levelname)s : %(message)s', datefmt='%m/%d/%Y %H:%M:%S')

bamfile, chromosomes_info = input_and_index(bamPath)
print("Counting unique reads")
total_unique_reads_count = count_unique_reads(bamfile, chromosomes_info)

# Effective genome length
L = args.e * sum(int(row[1]) for row in chromosomes_info)
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

f = open(outfile, 'wb')
for island in island_list:
    islandString = str(island[0]) + "\t" + str(island[1]) + "\t" + str(island[2]) + "\n"
    f.write(islandString)
# print(island_list[11104])
f.close()
bamfile.close()

print("Finished. Elapsed time, minutes: " + str((time.time() - startTime) / 60))
