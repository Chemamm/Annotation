#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--list_of_files", help="list of files from CVenrichment.py")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

Annotation.gene_overlapping(args.list_of_files, args.output)
