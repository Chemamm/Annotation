#!/usr/bin/env/python3

import argparse
import Annotation
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="CV file")
parser.add_argument("-o", "--output_label", help="output label")
parser.add_argument("-r", "--reference", default="", help="reference")
parser.add_argument("-p", "--pwd", help="path to working directory")
parser.add_argument("-t", "--target", help="target class to study")
parser.add_argument("-fo", "--final_output_label", help="output label for the final result")
parser.add_argument("-n", "--number", type=int, help="number of IDs to get")
parser.add_argument("-a", "--annotation_required", default=False, action="store_true", help="is the annotation required?")
parser.add_argument("-e", "--excel", default=False, help="complete excel")
parser.add_argument("-ef", "--expression_filter", default=False, type=float, help="FPKM to filter per condition")
args = parser.parse_args()


os.chdir(args.pwd)
os.system("mkdir %s" % args.target)
cvfile = Annotation.ConditionCV(args.file)
cvfile.getSpecificEnrichmentInputAnnotated(args.output_label, args.number, args.reference, args.excel, args.pwd, args.target, args.final_output_label, args.expression_filter, args.annotation_required)
