#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="trnscanse file")
parser.add_argument("-o", "--output", help="output filename label")
args = parser.parse_args()

Annotation.get_trna_gff(args.file, args.output)


