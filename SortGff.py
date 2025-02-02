#!/usr/bin/env/python3

from collections import defaultdict
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="fasta file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()


os.system("sort -k1,1 -k4,4n %s > tmp.gff" %args.file)
mRNA_exons = defaultdict(list)

# read the sorted gff file
with open('tmp.gff', 'r') as f:
    lines = f.readlines()

# group exons by mRNA parent
exon_dict = {}
for line in lines:
    fields = line.strip().split('\t')
    if fields[2] == 'mRNA':
        parent = fields[8].split(';')[0].split('=')[1]
        exon_dict[parent] = []
    elif fields[2] == 'exon':
        parent = fields[8].split(';')[1].split('=')[1]
        exon_dict[parent].append(fields)

# sort exons within each mRNA parent group
for parent, exons in exon_dict.items():
    sorted_exons = sorted(exons, key=lambda x: int(x[3]))
    exon_dict[parent] = sorted_exons

# write the sorted gff file with parenting mRNA first
with open(args.output, 'w') as f:
    for line in lines:
        fields = line.strip().split('\t')
        if fields[2] == 'mRNA':
            parent = fields[8].split(';')[0].split('=')[1]
            f.write(line)
            for exon in exon_dict[parent]:
                f.write('\t'.join(exon) + '\n')

os.system("rm tmp.gff")
