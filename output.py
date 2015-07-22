#!/usr/bin/env/ python


def write_output(outfile, island_list):
    f = open(outfile, 'wb')
    for island in island_list:
        island_string = str(island[0]) + "\t" + str(island[1]) + "\t" + str(island[2]) + "\n"
        f.write(island_string)
    # print(island_list[11104])
    f.close()

