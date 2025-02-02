#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="excel file")
parser.add_argument("-o", "--output_label", help="output label")
args = parser.parse_args()

cexcel = Annotation.CompleteExcel(args.file)
cexcel.getBackgrounds(args.output_label)
