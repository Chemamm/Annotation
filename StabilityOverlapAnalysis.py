#!/usr/bin/env/python3

import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f1", "--file1", help="stability overlap file")
parser.add_argument("-f2", "--file2", default=False, help="stability overlap file")
parser.add_argument("-d", "--degree", default=False, help="list of degrees separated by commas")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


sup_file = Annotation.SuperExactTest(args.file1)

if args.degree:
    degree_list = args.degree.split(",")
else:
    degree_list = False

sup_file.stability_analysis(args.output, degree_list, args.file2)
