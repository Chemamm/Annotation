#!/usr/bin/env python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

gff = Annotation.Gff(args.file)
gff.add_exon(args.output)