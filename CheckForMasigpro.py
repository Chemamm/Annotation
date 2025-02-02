#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--ematrix", help="counts ematrix")
parser.add_argument("-o", "--output", help="output ematrix")

args = parser.parse_args()

ematrix = Annotation.Ematrix(args.ematrix)
ematrix.check_for_masigpro(args.output)
