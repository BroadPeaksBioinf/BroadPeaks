#!/usr/bin/env python

import scipy
import numpy
"""
window = [start, number_of_tags]
island = [current_chromosome_name, island_start, window_start + window_size, island_score, island_number_of_reads, island_length - island_number_of_gaps, island_number_of_gaps])
"""


def calculate_and_append_score_for_fdr(island_list_treatment, window_list_control_dict, lambdaa, window_size, normalization):
    for island in island_list_treatment:
        current_chr = island[0]
        island_start = island[1]
        island_end = island[2]
        island_reads_count = island[3]

        island_length = island_end - island_start
        windows_per_island = island_length / window_size

        lambda_coeff = window_size / island_length * float(normalization)
        island_lambda = island_reads_count * lambda_coeff

        for i in range(len(window_list_control_dict[current_chr]) - windows_per_island):
            window_list = window_list_control_dict[current_chr]
            # first window in island
            first_window = window_list[i][0]
            last_window = window_list[i+windows_per_island][0]

            if island_start > first_window:
                continue
            elif island_start == first_window and island_end == last_window:
                # number of control tags per island
                control_reads_count = 0
                j = 0
                while j <= windows_per_island:
                    number_of_reads_per_window = window_list[i+j][1]
                    control_reads_count += number_of_reads_per_window
                    j += 1
                control_lambda = control_reads_count * lambda_coeff
                lambdaa = max(control_lambda, lambdaa)
                break
            break

        p = scipy.stats.poisson.pmf(island_lambda, lambdaa)
        if p < 1e-320:
            score_for_fdr = 1000
        else:
            score_for_fdr = -numpy.log(p)

        island.append(score_for_fdr)


def make_scores_dict(island_list):
    # {score : [island_position_in_list]}
    scores_dict = {}
    for i in range(len(island_list)):
        island = island_list[i]
        score_for_fdr = island[7]

        island_position_in_list = i
        # chr_name = island[0]

        try:
            scores_dict[score_for_fdr].append(island_position_in_list)
        except KeyError:
            scores_dict[score_for_fdr] = []
            scores_dict[score_for_fdr].append(island_position_in_list)
    return scores_dict


def calculate_and_append_fdr(island_list_treatment, island_list_control):
    treatment_scores_dict = make_scores_dict(island_list_treatment)
    treatment_scores = treatment_scores_dict.keys()
    treatment_scores.sort(reverse=True)

    control_scores = []
    for island in island_list_control:
        score_for_fdr = island[7]
        control_scores.append(score_for_fdr)
    control_scores.sort(reverse=True)

    negative_islands = 0
    true_islands = 0

    for t_score in treatment_scores:
        true_islands += 1
        while t_score <= control_scores[negative_islands]:
            negative_islands += 1

        fdr = 100.0 * negative_islands / true_islands

        for position in treatment_scores_dict[t_score]:
            island_list_treatment[position][8] = fdr
