#!/usr/bin/env/python3

import argparse
import Annotation
import os


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="expression matrix")
parser.add_argument("-o", "--output", help="output filename")
parser.add_argument("-ol", "--out_label", help="filename labels")
parser.add_argument("-n", "--n_filter", type=float, help="minimum cv to be considered unstable")
parser.add_argument("-m", "--m_filter", type=int, help="number of unstable transcripts to get")
parser.add_argument("-r", "--randomization", default=False, action="store_true", help="do the randomization")
args = parser.parse_args()

Annotation.stability(args.file, args.output, args.out_label, args.n_filter, args.m_filter, args.randomization)