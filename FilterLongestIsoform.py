#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file")
parser.add_argument("-o", "--output", help="output label")
args = parser.parse_args()

ft = Annotation.Fasta(args.fasta)
ft.filterlongestisoform(args.output)
