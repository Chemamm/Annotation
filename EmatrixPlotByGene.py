#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="list file")
parser.add_argument("-o", "--output", help="output folder")
parser.add_argument("-l", "--label",type=str, help="output folder")
parser.add_argument("-e", "--error",default=False, help="error sem matrix")
args = parser.parse_args()


em = Annotation.Ematrix(args.file)
if args.error ==False:
    em.expression_plot_by_gene(args.output, args.label)
else:
    em.expression_plot_by_gene_error(args.error, args.output,  args.label)