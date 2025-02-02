#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="gff file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()


Annotation.correct_rfam_gff(args.file, args.output)
