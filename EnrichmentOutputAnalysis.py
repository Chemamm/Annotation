#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--list_of_files", help="list of enrichment files")
parser.add_argument("-o", "--output_label", help="output label")
args = parser.parse_args()

Annotation.enrichment_output_analysis(args.list_of_files, args.output_label)
