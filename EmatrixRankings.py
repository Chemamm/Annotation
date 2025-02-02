import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--ematrix", help="Expression matrix")
parser.add_argument("-o", "--out_label", help="output filename")
parser.add_argument("-n", "--number", help="number of top genes to be analyzed")
args = parser.parse_args()


ematrix = Annotation.Ematrix(args.ematrix)
ematrix.get_ranking(args.out_label, args.number)
