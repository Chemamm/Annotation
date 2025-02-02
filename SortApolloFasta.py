#!/usr/bin/env/python3

import Annotation
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()


fa = Annotation.Fasta(args.file)
fa.sort_apollo_fasta(args.output)
