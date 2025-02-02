#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="family file")
parser.add_argument("-fa", "--fasta", help="fasta file")
parser.add_argument("-d", "--domain", help="domain file")
parser.add_argument("-o", "--out_label", help="output label filename")
args = parser.parse_args()


Annotation.defensins(args.fasta, args.file, args.domain, args.out_label)

