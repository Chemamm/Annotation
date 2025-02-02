#!/usr/bin/env/python3

import argparse
import Annotation
import os


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="CV file")
parser.add_argument("-o", "--output_label", help="output label")
parser.add_argument("-r", "--reference_label", default="", help="reference label")
parser.add_argument("-p", "--pwd", help="path to working directory")
parser.add_argument("-n", "--number", type=int, help="number of IDs to get")
parser.add_argument("-c", "--comparison", default=False, action="store_true", help="comparison background")
parser.add_argument("-e", "--excel", default=False, help="complete excel")
parser.add_argument("-ef", "--expression_filter", default=False, type=float, help="FPKM to filter per condition")
args = parser.parse_args()

os.chdir(args.pwd)
cv = Annotation.ConditionCV(args.file)
if args.comparison != False:
    cv.getEnrichmentInputAnnotated(args.output_label, args.number, args.reference_label, args.excel, args.pwd, args.expression_filter)
else:
    list_of_files = cv.getEnrichmentInput(args.output_label, args.number)
    Annotation.enrichment(list_of_files, args.reference_label, args.pwd, args.comparison, args.expression_filter)

