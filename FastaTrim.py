#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-o", "--output", help="output filename label")
parser.add_argument("-s", "--sep", default=" ", help="separator of name")
args = parser.parse_args()

fa = Annotation.Fasta(args.file)
fa.trim(args.sep, args.output)

