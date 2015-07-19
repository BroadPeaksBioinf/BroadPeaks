#!/usr/bin/env/ python

import sys
import argparse
parser = argparse.ArgumentParser(description="Tool for ChiP-seq analysis to find broad peaks")
parser.add_argument('-c', metavar='control', help="Path to `control.bam` file. DEFAULT: no control file", type=str)
args = parser.parse_args(['-c', 'st'])
print(args.control)

"""
try:
    options = parser.parse_args()
except:
    parser.print_help()
    sys.exit(0)
"""

#if "-h" in sys.argv or "--help" in sys.argv:
    #print "This is help text!"


"""
Now you see help for the test.py
"""
lst = [['i1', 1], ['i2', 2]]
for count in lst:
    print (count)
