#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-o", "--out_label", help="output name")
args = parser.parse_args()

fa = Annotation.Fasta(args.file)
fa.len_distribution(args.out_label)
