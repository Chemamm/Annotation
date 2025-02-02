#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file")
parser.add_argument("-i", "--ixoridb", help="IxoriDB browse output file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


fa = Annotation.Fasta(args.fasta)
ixori = Annotation.IxoriDBOut(args.ixoridb)
fa.filter(ixori.ids, args.output)
