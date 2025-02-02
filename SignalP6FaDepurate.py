#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fa file")
parser.add_argument("-s", "--signalp", help="signalp6 fasta file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


Annotation.signalp6_fa_depurate(args.signalp, args.file, args.output)
