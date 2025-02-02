#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cexcel", help="complete excel annotation file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

cexcel = Annotation.CompleteExcel(args.cexcel)
cexcel.get_annot_file(args.output)
