#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="expression matrix")
parser.add_argument("-o", "--out_label", help="output filename label")
parser.add_argument("-s", "--sep", help="separator of the sample name")
parser.add_argument("-m", "--method", default="ward", help="method for the hierarchical clustering")
parser.add_argument("-cn", "--color_number", type=int,  help="place of the string used for the color vector")
parser.add_argument("-fn", "--form_number", type=int, help="place of the string used for the form vector")
args = parser.parse_args()

ematrix = Annotation.Ematrix(args.file)
ematrix.hierarchical_clustering(args.out_label, args.method)
ematrix.correlation_matrix(args.out_label)
color_vector, form_vector = ematrix.get_vectors(args.sep, args.color_number, args.form_number)
ematrix.PCA(args.out_label, color_vector, form_vector)
ematrix.MDS(args.out_label, color_vector, form_vector)

