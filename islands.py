#!/usr/bin/env/ python

import numpy
import scipy


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
