#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="input label for files")
parser.add_argument("-fa", "--fasta", help="fasta file")
parser.add_argument("-g", "--gff", help="gff file")
parser.add_argument("-o", "--output", help="output label filename ")
args = parser.parse_args()


Annotation.getncRNA_gffcompare(args.file, args.output, args.gff, args.fasta)
