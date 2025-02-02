#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="expression matrix file")
parser.add_argument("-r", "--reference", help="reference genome file")
parser.add_argument("-o", "--output", help="output filename ")
args = parser.parse_args()



bg = Annotation.BlatGff(args.file)
bg.depurate_gff(args.reference, args.output)