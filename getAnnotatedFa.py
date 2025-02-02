#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-c", "--cexcel", help="complete excel")
parser.add_argument("-o", "--output", help="output name")
args = parser.parse_args()

cexcel = Annotation.CompleteExcel(args.cexcel)
fa = Annotation.Fasta(args.file)
fa.filter(cexcel.annotated, args.output)