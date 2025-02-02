import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Expression matrix")
parser.add_argument("-t", "--threshold", help="Threshold")
parser.add_argument("-s", "--samplesheet", help="Samplesheet file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

ematrix = Annotation.Ematrix(args.file)
ematrix.filterbyexpression_sampleSheet(args.threshold, args.samplesheet, args.output)
