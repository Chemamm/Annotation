#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-b", "--b2go", help="b2go table")
parser.add_argument("-e", "--excel", help="complete excel annotation file")
parser.add_argument("-o", "--obo", help="obo file")
args = parser.parse_args()

Annotation.numberofB2GO(args.excel, args.b2go, args.obo)
