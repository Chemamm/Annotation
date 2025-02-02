#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="ixoridb file")
parser.add_argument("-b", "--bed", help="bed file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

bed = Annotation.Bed(args.bed)
ixoridb = Annotation.IxoriDBOut(args.file)
bed.intersect(ixoridb.ids, args.output)