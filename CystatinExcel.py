#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-o", "--output", help="output filename label")
args = parser.parse_args()

Annotation.excel_cystatin(args.file, args.output)


