#!/usr/bin/env/ python

import logging
import os.path
import pysam


def get_chromosomes_info(bam_path):
    # Check if there is an index file, create one if there isn't
    if not os.path.isfile(bam_path + ".bai"):
        pysam.index(bam_path)
        logging.info('No BAM index file was found, new index was generated : `{}`'.format(bam_path + ".bai"))
    # Take chromosome data from BAM index:
    # (ref.seq. name, ref.seq. length, number of mapped reads and number of unmapped reads)
    chromosomes_info = []
    logging.info('Collecting information about sample from .bai file: '
                 '[ref.seq. name, ref.seq. length, number of mapped and unmapped reads]')
    for chr in pysam.idxstats(bam_path):
        chromosomes_info.append(chr.split("\t")[:-1])
    # Last line is unmapped reads, we don't need them
    chromosomes_info.pop()
    # print(chromosomes_info)
    return chromosomes_info


def count_unique_reads(bam_path, chromosomes_info):
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    total_unique_reads_count = 0
    previous_read_strand = 0
    for chromosome in chromosomes_info:
        chr_unique_reads_count = 0
        chr_total_reads_count = 0
        beginning_of_the_previous_read = 0
        current_chromosome_name = chromosome[0]
        logging.info("Counting reads on {}".format(current_chromosome_name))
        # currentChromosomeSize = int(chromosome[1])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
        for read in all_reads_in_chromosome:
            read_str = str(read)
            # read strand: 0 = +         16 = -
            read_strand = ([int(s) for s in read_str.split() if s.isdigit()][0])
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            if beginning_of_the_read != beginning_of_the_previous_read or (read_strand != previous_read_strand):
                beginning_of_the_previous_read = beginning_of_the_read
                previous_read_strand = read_strand
                total_unique_reads_count += 1
                chr_unique_reads_count += 1
            chr_total_reads_count += 1
        logging.info("There are {} unique reads among {} on {}".format(chr_unique_reads_count,
                                                                       chr_total_reads_count,
                                                                       current_chromosome_name))
    # print("Unique reads counted")
    bamfile.close()
    return total_unique_reads_count
    # normalising_coefficient = total_unique_reads_count / 1000000
    # it can help to calculate experiments with "control"
    # read_coverage has to be multiplied on normalising_coefficient


def count_effective_length(effective_proportion, chromosomes_info):
    total_genome_length = sum(int(row[1]) for row in chromosomes_info)
    effective_length = effective_proportion * total_genome_length
    return effective_length


def count_lambda(unique_reads_count, window_size, effective_length):
    lambdaa = float(window_size) * float(unique_reads_count) / float(effective_length)
    return lambdaa
