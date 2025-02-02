#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Spearman correlation results file")
parser.add_argument("-c", "--cexcel", help="Complete Excel")
parser.add_argument("-o", "--output_label", help="Output filename label")
parser.add_argument("-t", "--tissue", help="Tissue")
parser.add_argument("-pl", "--plot_label", help="Label for the distribution plot")
parser.add_argument("-nb", "--nbins", type=int, default=20, help="Number of bins to plot")
parser.add_argument("-n", "--number", type=int, default=200, help="Number of genes to get as top or bot")
args = parser.parse_args()

Annotation.spearman_corr_byweight_distribution(args.file, args.output_label, args.plot_label, args.nbins)
Annotation.spearman_corr_byweight_top(args.file, args.cexcel, args.tissue, args.output_label, args.number)