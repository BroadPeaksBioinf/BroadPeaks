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
import arguments
import pre_counting
import islands
import output


startTime = time.time()

# ALL TO main()

parser = argparse.ArgumentParser(description="Tool for ChiP-seq analysis to find broad peaks")

parser.add_argument('infile', help="Path to `input.bam` file", type=str)
parser.add_argument('-w', dest='window_size', help="Window size (bp).  DEFAULT: 200", type=int, default=200)
parser.add_argument('-g', dest='gap', help="Gap size shows how many windows could be skipped. DEFAULT: 1",
                    type=int, default=1, choices=[1, 2, 3])
parser.add_argument('-t', dest='threshold', help="Island score threshold. DEFAULT: 100", type=int, default=0)
parser.add_argument('-o', dest='outdir', help="Path to directory for output `*_peaks.bed` file. "
                                              "DEFAULT: output will be in the same directory as `input.bam`",
                    type=str)
parser.add_argument('-n', dest="output_name", help="Specify output name. "
                                                   "DEFAULT : an input file name + `_peaks.bed`", type=str)
parser.add_argument('-e', help="Proportion of effective genome length; has to be in range(0.0, 1.0) DEFAULT: 0.77",
                    type=float, default=0.77)
parser.add_argument('-c', dest='control', help="Path to `control.bam` file. DEFAULT: no control file",
                    type=str)
parser.add_argument('--log', help="To see only current run LOG file. "
                                  "DEFAULT : LOG file contains information from all runs", action='store_true')
# parser.add_argument('-p', dest='p_value', help="p-value; has to be in range(0.0, 1.0). DEFAULT: 0.01", type=float, default=0.01)
# parser.add_argument('-ref', dest='reference genome', help="Reference genome.  DEFAULT: 'hg19'", type=str, default='hg19')

"""
INPUT ARGUMENTS:

see broadPeaks.py -h or --help
bamPath = "/home/dima/BAMfiles/Bernstein_H1_hESC_CTCF.bam"
"/home/dima/BAMfiles/h3k4me3_rep1.bam"
'/home/yegor/Alex_project/H3K4me3.bam'
"""
# args as list of strings
#args = parser.parse_args(['/home/yegor/Alex_project/H3K4me3.bam'])
# ["/home/user/SICERproj/BAMfiles/H3K4Me3_test.bam"
args = parser.parse_args()

bamPath = arguments.check_input(args.infile)
arguments.make_log(bamPath, args.log)
WINDOW_SIZE = args.window_size
GAP = args.gap
EFFECTIVE_PROPORTION = arguments.check_effective_proportion(args.e)
ISLAND_SCORE_THRESHOLD = args.threshold
outfile = arguments.check_outfile(args.outdir, args.output_name, bamPath)
# controlPath = arguments.check_control(args.control)
# p0 = arguments.check_p_value(args.p_value)
p0 = 0.01

# main_functions
chromosomes_info = pre_counting.get_chromosomes_info(bamPath)

logging.info("\nStep 1 of 4\nCOUNTING UNIQUE READS\n")
total_unique_reads_count = pre_counting.count_unique_reads(bamPath, chromosomes_info)

# Effective genome length (L)
effective_length = pre_counting.count_effective_length(EFFECTIVE_PROPORTION, chromosomes_info)
# Lambda for poisson distribution
lambdaa = pre_counting.count_lambda(total_unique_reads_count, WINDOW_SIZE, effective_length)
# Minimum #reads in a window for eligibility
# Formula (1), finding l0
# Must make more clear variable name
l0 = scipy.stats.poisson.ppf(1 - p0, lambdaa)
logging.info("\nWindow read threshold is {} reads, \ni.e. {} is minimum number of reads in window "
             "to consider this window `eligible` with Poisson distribution p-value {}".format(l0, l0, p0))

logging.info("\nStep 2 of 4\nMAKING WINDOW LIST\n")
window_list = islands.make_windows_list(bamPath, chromosomes_info, l0, WINDOW_SIZE, GAP, total_unique_reads_count)

logging.info("\nStep 3 of 4\nMAKING ISLAND LIST\n")
island_list = islands.make_islands_list(window_list, lambdaa, WINDOW_SIZE, l0, chromosomes_info, ISLAND_SCORE_THRESHOLD)

logging.info("\nStep 4 of 4\nWRITING FOUND ISLANDS TO `{}` BED FILE\n".format(outfile))
output.write_output(outfile, island_list)
# print(len(windowList), sys.getsizeof(windowList)/1024)

"""
IF WE HAVE CONTROL

control_bam = check_and_open_input_bam(controlPath, LOG_filename)
control_chromosomes_info = get_chromosomes_info(controlPath)
control_unique_reads_count = count_unique_reads(control_bam, control_chromosomes_info)
control_L = count_effective_length(EFFECTIVE_PROPORTION, control_chromosomes_info)
control_lambda = count_lambda(control_unique_reads_count)
control_l0 = scipy.stats.poisson.ppf(1 - p0, control_lambda)
control_windowList = make_windows_list(control_bam, control_chromosomes_info, control_l0, WINDOW_SIZE, GAP,
                                       control_unique_reads_count)
"""

print("Finished. Elapsed time, minutes: " + str((time.time() - startTime) / 60))
