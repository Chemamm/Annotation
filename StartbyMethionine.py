import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

Annotation.startbymethionine(args.fasta, args.output)
