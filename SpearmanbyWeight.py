#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="complete excel file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

Annotation.spearman_corr_orderedbyweight(args.file, args.output)
