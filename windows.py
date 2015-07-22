#!/usr/bin/env/ python
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
                            gap_count = 0
                        # * unique_reads_count/1000000 is for normalization per million reads
                        # now we are able to compare control and sample

                        # NEED TO CHANGE LAMBDA AND N
                        window_list.append([i, window_reads_count])
                        # / (unique_reads_count/1000000)])
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
    bamfile.close()
    return window_list