#!/usr/bin/env/ python
import logging

def write_output(outfile, island_list, ISLAND_SCORE_THRESHOLD):
    f = open(outfile, 'wb')
    logging.info("Line in BED file: 'chromosome_name\tisland_start\tisland_end\tisland_score\tnumber_of_reads_in_the_island\tisland_length_in_windows\tnumber_of_gaps_in_the_island'")
    for island in island_list:
        if island[3] >= ISLAND_SCORE_THRESHOLD:
            island_string = str(island[0]) + "\t" + str(island[1]) + "\t" + str(island[2]) + "\t" + str(island[3]) + "\t" + str(island[4])+ "\t" + str(island[5])+ "\t" + str(island[6]) + "\n"
            f.write(island_string)
    # print(island_list[11104])
    f.close()
