#!/usr/bin/env python3

import argparse, sys
# import math, time, re
import gzip
import numpy as np
import pandas as pd
# from scipy import stats
# from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (cchiang@genome.wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-03-27 09:43 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
extract_LD_class.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: extract plink LD values based on variant attributes")
    parser.add_argument('-i', '--input',
                        metavar='FILE', dest='input_path',
                        type=str, default=None,
                        help='plink file of LD [stdin]')
    parser.add_argument('-a', '--attrib',
                        metavar='FILE', dest='attrib_path',
                        required=True,
                        type=str, default=None,
                        help='tab-delimited file of variant attributes')
    parser.add_argument('--cadd_a',
                        metavar='STR', dest='cadd_a',
                        required=False,
                        type=str, default='0,1000',
                        help='CADD phred score range A. E.g. "0,20". [0,1000]')
    parser.add_argument('--cadd_b',
                        metavar='STR', dest='cadd_b',
                        required=False,
                        type=str, default='0,1000',
                        help='CADD phred score range B. E.g. "0,20".[0,1000]')

    # parse the arguments
    args = parser.parse_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)

    # send back the user input
    return args

# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data    

def scalarize(x):
    if isinstance(x, pd.Series):
        return(x.iloc[0])
    else:
        return(x)

def listize(x):
    if isinstance(x, pd.Series):
        return(x)
    else:
        return [x]

# check if value (float) is in range of 2-tuple [float, float]
def in_range(x, myrange):
    return (x >= myrange[0] and x < myrange[1])    

# primary function
def myfunc(input_file,
           attrib_path,
           cadd_range_a,
           cadd_range_b):
    vardf = pd.read_csv(attrib_path, sep='\t',
                        index_col=2, # variant id is row index
                        names=['chr', 'pos', 'id', 'gene', 'func', 'caddphred'])
    plink_header = input_file.readline().rstrip().split()
    print('\t'.join(plink_header + ['CADD_A', 'CADD_B']))

    # i = 0
    for line in input_file:
        # i = i + 1
        # if (not i % 10000):
        #     sys.stderr.write("%s\n" % i)
        v = line.rstrip().split()
        varpair = dict(zip(plink_header, v))

        # -----------------------------------
        # at beginning, set pass flag to 0 (false)
        var_pass = 0

        # -----------------------------------
        # check if variants are in the same gene
        gene_a = listize(vardf.loc[varpair['SNP_A'], 'gene'])
        gene_b = listize(vardf.loc[varpair['SNP_B'], 'gene'])

        # print('a', set(gene_a), 'b', set(gene_b))

        if (set(gene_a).isdisjoint(set(gene_b))):
            # print('disjoint!')
            continue # go to next iteration

        # -----------------------------------
        # check if CADD scores are within range
        cadd_a = scalarize(vardf.loc[varpair['SNP_A'], 'caddphred'])
        cadd_b = scalarize(vardf.loc[varpair['SNP_B'], 'caddphred'])

        if in_range(cadd_a, cadd_range_a) and in_range(cadd_b, cadd_range_b):
            var_pass = 1
        elif in_range(cadd_a, cadd_range_b) and in_range(cadd_b, cadd_range_a):
            var_pass = 1

        if var_pass:
            print('\t'.join(map(str, v + [cadd_a, cadd_b])))

    # print(vardf)
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.input_path == None:
        input_file = sys.stdin
    else:
        input_file = get_file(args.input_path)

    # # get permutation data
    # permutation_file = get_file(args.permutation_path)

    cadd_range_a = list(map(float, args.cadd_a.strip().split(',')))
    cadd_range_b = list(map(float, args.cadd_b.strip().split(',')))

    # call primary function
    myfunc(input_file,
           args.attrib_path,
           cadd_range_a,
           cadd_range_b)

    # close the files
    input_file.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
