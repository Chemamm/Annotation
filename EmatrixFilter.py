#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Expression matrix")
parser.add_argument("-c", "--cexcel", default=False, help="complete excel annotation file")
parser.add_argument("-o", "--output", help="output filename")
parser.add_argument("-f2", "--file2", default=False, help="Expression matrix FPKM")
parser.add_argument("-t", "--threshold", default=False, type=int, help="FPKM threshold")
parser.add_argument("-n", "--nconditions", default=False, type=int, help="Number of conditions")
parser.add_argument("-r", "--replicates", default=False, type=int, help="number of replicates")
args = parser.parse_args()

if args.file2 == False and args.cexcel != False:
    em = Annotation.Ematrix(args.file)
    cexcel = Annotation.CompleteExcel(args.cexcel)
    em.filter(cexcel.transcriptIDs, args.output)
else:
    em1 = Annotation.Ematrix(args.file)
    em2 = Annotation.Ematrix(args.file2)
    IDs = em2.filterbyexpression(args.replicates, args.nconditions, args.threshold)
    em1.filter(IDs, args.output)
