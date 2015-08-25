#!/usr/bin/env python

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


def broadpeaks_with_control(bam_path, control_path, window_size, gap, EFFECTIVE_PROPORTION, ISLAND_SCORE_THRESHOLD, p0):

    chromosomes_info = pre_counting.get_chromosomes_info(bam_path)
    control_chromosomes_info = pre_counting.get_chromosomes_info(control_path)

    logging.info("\nStep 1 of 4\nCOUNTING UNIQUE READS\n")
    input_unique_reads_count = pre_counting.count_unique_reads(bam_path, chromosomes_info)
    control_unique_reads_count = pre_counting.count_unique_reads(control_path, control_chromosomes_info)

    # Effective genome length (L)
    effective_length = pre_counting.count_effective_length(EFFECTIVE_PROPORTION, chromosomes_info)

    # Lambda for poisson distribution
    lambdaa = pre_counting.count_lambda(input_unique_reads_count, window_size, effective_length)

    # Minimum #reads in a window for eligibility
    # Formula (1), finding l0
    l0 = scipy.stats.poisson.ppf(1 - p0, lambdaa)
    logging.info("\nWindow read threshold is {} reads, \ni.e. {} is minimum number of reads in window "
                 "to consider this window `eligible` with Poisson distribution p-value {}".format(l0, l0, p0))


    logging.info("\nStep 2 of 4\nMAKING WINDOW LIST\n")
    NORMALIZATION_CONSTANT = float(input_unique_reads_count)/float(control_unique_reads_count)
    window_list_input = islands.make_windows_list(bam_path, chromosomes_info, l0, window_size, gap, input_unique_reads_count, NORMALIZATION_CONSTANT)
    window_list_control = islands.make_windows_list(bam_path, chromosomes_info, l0, window_size, gap, control_unique_reads_count, NORMALIZATION_CONSTANT)

    #window_list = islands.modify_window_list_based_on_control(control_path, chromosomes_info, l0, window_size, gap, input_unique_reads_count, control_unique_reads_count, window_list_temp)


    logging.info("\nStep 3 of 4\nMAKING ISLAND LIST\n")
    island_list_input = islands.make_islands_list(window_list_input, lambdaa, window_size, l0, chromosomes_info, ISLAND_SCORE_THRESHOLD)
    island_list_control = islands.make_islands_list(window_list_control, lambdaa, window_size, l0, chromosomes_info, ISLAND_SCORE_THRESHOLD)

    island_list = islands.find_unintersected_islands(island_list_input,island_list_control)

    # calculate FDR
    FDR = (len(island_list_control) - (len(island_list_input)-len(island_list)))/len(island_list)
    logging.info("\nFDR is {} reads, \n".format(FDR))

    return(island_list)