#!/usr/bin/env/ python

import sys
import logging
import os.path
import colorized_log


def check_input(input_path):
    """
    :param input_path: (str)
    :return: input_path
    """
    # Check if there is a file at entered dir
    if not os.path.isfile(input_path):
        sys.exit("No BAM file specified in input '{}' or there is no such a file".format(input_path))
    # Check if it is a BAM file
    if input_path[-4:] != '.bam':
        sys.exit("`{}` is other file type, not BAM. This tool works only with BAM-files as input (*.bam)".
                 format(input_path))
    return input_path


def make_log(input_path, log_flag):
    # log_path = input_path[:-4]
    log_path = os.path.dirname(input_path) + "/" + os.path.basename(input_path)[:-4] + "_log.log"
    if log_flag:
        logging.basicConfig(filename=log_path, level=logging.DEBUG,
                            format='%(asctime)s : %(levelname)s : %(message)s', datefmt='%m/%d/%Y %H:%M:%S',
                            filemode='w')
    else:
        logging.basicConfig(filename=log_path, level=logging.DEBUG,
                            format='%(asctime)s : %(levelname)s : %(message)s', datefmt='%m/%d/%Y %H:%M:%S')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # console.setFormatter(formatter)
    formatter = colorized_log.ColoredFormatter('%(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


def check_control(control_path):
    if control_path != "unspecified":
        if not os.path.isfile(control_path):
            logging.error("No control BAM file specified in input '{}' or there is no such a file".format(control_path))
            sys.exit("`{}` is not a path to BAM file. \n More information in LOG file for this run".format(control_path))
        if control_path[-4:] != '.bam':
            logging.error("`{}` is other file type, not BAM. This tool works only with BAM-files as control input (*.bam)".
                          format(control_path))
            sys.exit("`{}` is not a BAM file. \n More information in `{}`".format(control_path))
    return control_path


def check_p_value(p0):
    if p0 <= 0 or p0 >= 1:
        logging.error("`{}` is incorrect p-value. p-value has to be in range(0.0, 1.0)".format(p0))
        sys.exit()
    return p0


def check_effective_proportion(effective_proportion):
    if effective_proportion <= 0 or effective_proportion > 1:
        logging.error("`{}` is incorrect proportion of effective genome length. \n "
                      "Proportion of effective genome length has to be in range(0.0, 1.0)".format(effective_proportion))
        sys.exit()
    return effective_proportion


def check_outfile(output_dir, output_name, input_path):
    if not output_dir:
        output_dir = os.path.dirname(input_path)
    if not output_name:
         output_name = os.path.basename(input_path)[:-4] + "_peaks"
    # must test validity of output_name as filename
    outfile = output_dir + "/" + output_name + ".bed"
    return outfile
