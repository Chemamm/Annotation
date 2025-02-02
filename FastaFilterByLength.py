import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Apriori intersecter file")
parser.add_argument("-min", "--min", default=0, type=float, help="Minimum length")
parser.add_argument("-max", "--max", default=100000000000000, type=float, help="Maximum length")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

fa = Annotation.Fasta(args.file)
fa.filter_by_length(args.output, args.max, args.min)