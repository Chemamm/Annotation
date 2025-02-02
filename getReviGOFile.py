#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="enrichment output matrix")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

en_file = Annotation.EnrichmentOut(args.file)
en_file.getReviGOFile(args.output)