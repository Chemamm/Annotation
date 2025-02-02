#!/usr/bin/env/python3

import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_label", help="output label name")
args = parser.parse_args()

Annotation.weight_figure(args.output_label)
