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
                        print ('again wrong')
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
