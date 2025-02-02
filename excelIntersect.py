#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file")
parser.add_argument("-e", "--excel", help="complete excel annotation file")
parser.add_argument("-o", "--output", help="output filename")
parser.add_argument("-r", "--reverse", action="store_true", default=False, help="reverse intersection")
args = parser.parse_args()

Annotation.intersectexcel(args.excel, args.fasta, args.output, args.reverse)
