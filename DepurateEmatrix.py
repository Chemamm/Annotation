#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="expression matrix file")
parser.add_argument("-fi", "--filter", default=False, type=int, help="minimum TPM/FPKM")
parser.add_argument("-n", "--nreplicates", default=False, help="number of replicates")
parser.add_argument("-m", "--matrix", default=False, help="expression matrix file for filtering")
args = parser.parse_args()


em = Annotation.Ematrix(args.file)
em.depurate(args.filter, args.nreplicates, args.matrix)
