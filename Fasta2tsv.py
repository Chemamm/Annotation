#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-fa", "--fasta", help="fasta file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

fa = Annotation.Fasta(args.fasta)
fa.fasta2tsv(args.output)

