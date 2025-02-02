#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="gff file")
parser.add_argument("-fa", "--fasta", help="fasta file")
parser.add_argument("-l", "--label", help="gff file")
parser.add_argument("-o", "--output", help="output filename ")
args = parser.parse_args()

gff = Annotation.Gff(args.file)
gff.get_sequences_from_IDs(args.label, args.fasta, args.output)

