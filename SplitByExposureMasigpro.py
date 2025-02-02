#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--ematrix", help="normalized counts ematrix with unfed")
parser.add_argument("-o", "--output_label", help="output label")
parser.add_argument("-t", "--tissue", help="MG or SG")
parser.add_argument("-ku", "--keep_unfed", default=False, action="store_true", help="keeping unfeds")
args = parser.parse_args()

ematrix = Annotation.Ematrix(args.ematrix)
ematrix.split_by_exposure_masigpro(args.output_label, args.tissue, args.keep_unfed)