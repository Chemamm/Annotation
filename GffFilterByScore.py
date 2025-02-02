#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="gff file")
parser.add_argument("-s", "--score", type=float, help="minimum score")
parser.add_argument("-o", "--output", help="output filename")

args = parser.parse_args()

gff = Annotation.Gff(args.file)
gff.filter_by_score(float(args.score), args.output)

