#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="input label for files")
parser.add_argument("-t", "--type", help="class for tsv")
parser.add_argument("-o", "--output", help="output label filename ")
args = parser.parse_args()

gff = Annotation.Gff(args.file)
gff.to_tsv(args.type, args.output)

