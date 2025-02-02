#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="list file")
parser.add_argument("-fa", "--fasta", help="fasta file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


with open(args.file) as input:
    input = list(input)

filter_list = []
for line in input:
    filter_list.append(line.replace("\n",""))

fa = Annotation.Fasta(args.fasta)
fa.filter(filter_list, args.output)