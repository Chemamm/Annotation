#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--ematrix", help="counts ematrix")
parser.add_argument("-o", "--output_label", help="output label filename")
args = parser.parse_args()

ematrix = Annotation.Ematrix(args.ematrix)
ematrix.rename_and_split_masigpro(args.output_label)