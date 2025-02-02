#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-e", "--ematrix", help="expression matrix")
parser.add_argument("-o", "--out_label", help="output filename label")
parser.add_argument("-t", "--threshold", type=float, help="expression threshold")
args = parser.parse_args()



em = Annotation.Ematrix(args.ematrix)
em.filterbyexpression_MGSG(args.file, args.threshold, args.out_label)

