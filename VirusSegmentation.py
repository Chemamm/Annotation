#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="genome fasta file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

Annotation.virus_segmentation(args.file, args.output)

