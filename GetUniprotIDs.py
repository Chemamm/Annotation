import Annotation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="DE file")
parser.add_argument("-e", "--excel", help="Complete MGSG Annotation tsv/excel")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

gene_list = []
with open(args.file) as genes:
    lines = list(genes)

for line in lines:
    gene = line.replace("\n","")
    if gene != "" and gene != "x":
        gene_list.append(gene)

annot = Annotation.CompleteExcel(args.excel)
annot.getUniprotIDs(gene_list, args.output)