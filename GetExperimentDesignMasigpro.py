#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--ematrix", help="counts ematrix")
parser.add_argument("-o", "--output_label", help="output label filename")
parser.add_argument("-r", "--replicates", help="number of replciates")

args = parser.parse_args()

ematrix = Annotation.Ematrix(args.ematrix)
ematrix.get_experiment_design_masigpro(args.output_label, args.replicates)