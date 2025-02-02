#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="expression matrix")
parser.add_argument("-o", "--output", help="output filename label")
parser.add_argument("-c", "--cexcel", help="Complete Excel filename")
args = parser.parse_args()

ematrix = Annotation.Ematrix(args.file)
ematrix.genestoqPCR(args.cexcel, args.output)
