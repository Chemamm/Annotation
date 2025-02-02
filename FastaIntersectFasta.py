#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file")
parser.add_argument("-f2", "--fasta2", help="fasta file")
parser.add_argument("-o", "--output", help="output filename")
parser.add_argument("-r", "--reverse", action="store_true", default=False, help="reverse intersection")
args = parser.parse_args()

fa = Annotation.Fasta(args.fasta)
fa.intersect(args.fasta2, args.output, args.reverse)
