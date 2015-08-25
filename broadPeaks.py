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
import broadpeaks_wo_control
import broadpeaks_with_control

startTime = time.time()

# ALL TO main()

parser = argparse.ArgumentParser(description="Tool for ChiP-seq analysis to find broad peaks")

parser.add_argument('infile', help="Path to `input.bam` file", type=str)
parser.add_argument('-w', dest='window_size', help="Window size (bp).  DEFAULT: 200", type=int, default=200)
parser.add_argument('-g', dest='gap', help="Gap size shows how many windows could be skipped. DEFAULT: 1",
                    type=int, default=1, choices=[1, 2, 3])
parser.add_argument('-t', dest='threshold', help="Island score threshold. DEFAULT: 0", type=int, default=0)
parser.add_argument('-o', dest='outdir', help="Path to directory for output `*_peaks.bed` file. "
                                              "DEFAULT: output will be in the same directory as `input.bam`",
                    type=str)
parser.add_argument('-n', dest="output_name", help="Specify output name. "
                                                   "DEFAULT : an input file name + `_peaks.bed`", type=str)
parser.add_argument('-e', help="Proportion of effective genome length; has to be in range(0.0, 1.0) DEFAULT: 0.77",
                    type=float, default=0.77)
parser.add_argument('-c', dest='control', help="Path to `control.bam` file. DEFAULT: no control file",
                    type=str, default="unspecified")
parser.add_argument('--log', help="To see only current run LOG file. "
                                  "DEFAULT : LOG file contains information from all runs", action='store_true')
# parser.add_argument('-p', dest='p_value', help="p-value; has to be in range(0.0, 1.0). DEFAULT: 0.01", type=float, default=0.01)
# parser.add_argument('-ref', dest='reference genome', help="Reference genome.  DEFAULT: 'hg19'", type=str, default='hg19')

"""
INPUT ARGUMENTS:

see broadPeaks.py -h or --help
bam_path = "/home/dima/BAMfiles/Bernstein_H1_hESC_CTCF.bam"
"/home/dima/BAMfiles/h3k4me3_rep1.bam"
'/home/yegor/Alex_project/H3K4me3.bam'
"""
# args as list of strings
# args = parser.parse_args(['/media/user/DISK1/SICER_project/Inputs_mouse/GSM1288312.bam'])
# ["/home/user/SICERproj/BAMfiles/H3K4Me3_test.bam"
args = parser.parse_args(['/media/user/DISK1/SICER_project/Inputs_mouse/GSM1562338.bam'])

bam_path = arguments.check_input(args.infile)
arguments.make_log(bam_path, args.log)
WINDOW_SIZE = args.window_size
GAP = args.gap
EFFECTIVE_PROPORTION = arguments.check_effective_proportion(args.e)
ISLAND_SCORE_THRESHOLD = args.threshold
outfile = arguments.check_outfile(args.outdir, args.output_name, bam_path)
control_path = arguments.check_control(args.control)
# p0 = arguments.check_p_value(args.p_value)
p0 = 0.1

# main_functions
if control_path == "unspecified":
    island_list = broadpeaks_wo_control.broadpeaks_wo_control(bam_path, WINDOW_SIZE, GAP,
                                                              EFFECTIVE_PROPORTION, ISLAND_SCORE_THRESHOLD, p0)
else:
    island_list = broadpeaks_with_control.broadpeaks_with_control(bam_path, control_path, WINDOW_SIZE, GAP,
                                                                  EFFECTIVE_PROPORTION, ISLAND_SCORE_THRESHOLD, p0)


logging.info("\nStep 4 of 4\nWRITING FOUND ISLANDS TO `{}` BED FILE\n".format(outfile))
output.write_output(outfile, island_list, ISLAND_SCORE_THRESHOLD)


print("Finished. Elapsed time, minutes: " + str((time.time() - startTime) / 60))
