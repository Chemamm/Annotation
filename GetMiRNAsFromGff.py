#!/usr/bin/env/python3

import Annotation
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="gff file")
parser.add_argument("-o", "--out_label", help="output file label")
args = parser.parse_args()


Annotation.get_miRNAs_from_gff(args.file, args.out_label)
