
import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Ematrix file")
parser.add_argument("-o", "--output", help="output")
args = parser.parse_args()

Annotation.get_MGSG_means(args.file, args.output)


