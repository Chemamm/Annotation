#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="complete excel file")
parser.add_argument("-o", "--output", help="output label")
args = parser.parse_args()

excel = Annotation.CompleteExcel(args.file)
excel.splitSecreted(args.output)
