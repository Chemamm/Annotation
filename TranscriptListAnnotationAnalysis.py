import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="transcript list file")
parser.add_argument("-c", "--cexcel", help="complete excel file")
parser.add_argument("-o", "--output_label", help="output label name")
parser.add_argument("-s", "--sep", default=False, help="separator of the list ")
parser.add_argument("-d", "--depurate", default=False, action="store_true", help="whether to split by . transcript names or not")
args = parser.parse_args()

tl = Annotation.TranscriptList(args.file, args.sep)
if args.depurate:
    tl.depurate()

tl.annotationAnalysis(args.cexcel, args.output_label)