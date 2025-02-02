import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--excel", help="complete excel")
parser.add_argument("-o", "--output", help="output label")
args = parser.parse_args()

ex = Annotation.CompleteExcel(args.excel)
ex.expression_distribution(args.output)
