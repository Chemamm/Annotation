import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="background file")
parser.add_argument("-ob", "--obo", help="GO obo filename")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

gob = Annotation.GOBackground(args.file)
gob.get_ancestors(args.obo, args.output)
