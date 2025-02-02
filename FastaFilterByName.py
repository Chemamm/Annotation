#!/usr/bin/env/python3

import Annotation
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-n", "--name", help="string to filter")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()


fa = Annotation.Fasta(args.file)
fa.filter_by_name(args.name, args.output)
