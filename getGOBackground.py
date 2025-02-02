import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cexcel", help="complete excel file")
parser.add_argument("-o", "--output", help="output name")
args = parser.parse_args()

cexcel = Annotation.CompleteExcel(args.cexcel)
cexcel.getGOBackground(args.output)
