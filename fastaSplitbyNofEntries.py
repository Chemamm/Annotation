#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file")
parser.add_argument("-n", "--n", help="number of entries the file must have")
parser.add_argument("-o", "--output", help="output label")
args = parser.parse_args()

ft = Annotation.Fasta(args.fasta)
ft.split(args.n, args.output)
