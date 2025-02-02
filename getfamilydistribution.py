import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cexcel", help="complete excel annotation file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

Annotation.family_distribution_stats(args.cexcel, args.output)
