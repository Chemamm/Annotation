#!/usr/bin/env/python3

import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="DE overlap file")
parser.add_argument("-o", "--output", help="output name")
args = parser.parse_args()

sup_file = Annotation.SuperExactTest(args.file)
sup_file.DE_analysis(args.output)
