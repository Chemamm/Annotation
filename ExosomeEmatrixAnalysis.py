#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="expression matrix file")
parser.add_argument("-o", "--out_label", help="output filename label")
args = parser.parse_args()



Annotation.exosome_ematrix_analysis(args.file, args.out_label)
