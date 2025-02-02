import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="fasta file")
args = parser.parse_args()


fa = Annotation.Fasta(args.fasta)
fa.print_lengths()