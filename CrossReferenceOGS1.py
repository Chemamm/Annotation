#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f1", "--file1", help="fa file")
parser.add_argument("-f2", "--file2", help="fa file")
parser.add_argument("-g", "--gff", help="annotation file")
parser.add_argument("-o", "--output", help="output filename ")
args = parser.parse_args()

Annotation.crossreference_OGS1(args.file1, args.file2, args.gff, args.output)
