#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-e", "--ematrix", help="expresion matrix")
parser.add_argument("-o", "--output", help="output name")
args = parser.parse_args()

fa = Annotation.Fasta(args.file)
ematrix = Annotation.Ematrix(args.ematrix)
ematrix.filter(fa.IDs, args.output)
