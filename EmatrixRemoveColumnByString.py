#!/usr/bin/env/python3

import Annotation
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="ematrix file")
parser.add_argument("-s", "--string", help="string to filter the df")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()


em = Annotation.Ematrix(args.file)
em.remove_columns_containing_string(args.string, args.output)


