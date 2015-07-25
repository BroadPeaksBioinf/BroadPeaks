#!/usr/bin/env/ python
import logging

def write_output(outfile, island_list):
    f = open(outfile, 'wb')
    logging.info("Line in BED file: 'chromosome_name\tisland_start\tisland_end\tisland_score'")
    for island in island_list:
        island_string = str(island[0]) + "\t" + str(island[1]) + "\t" + str(island[2]) + "\t" + str(island[3]) + "\n"
        f.write(island_string)
    # print(island_list[11104])
    f.close()
