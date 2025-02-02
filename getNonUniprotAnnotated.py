#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="excel best alignment file")
parser.add_argument("-u", "--uniprot", help="uniprot annotation file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

Annotation.getUnknown(args.file, args.uniprot, args.output)
