#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--ematrix", help="counts ematrix")
parser.add_argument("-o", "--output_label", help="output label filename")
parser.add_argument("-ku", "--keep_unfed", default=False, action="store_true", help="output ematrix")
args = parser.parse_args()

ematrix = Annotation.Ematrix(args.ematrix)
ematrix.split_by_tissue_masigpro(args.output_label, args.keep_unfed)