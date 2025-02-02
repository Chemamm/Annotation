#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="list of IDs")
parser.add_argument("-o", "--output", help="output filename label")
parser.add_argument("-c", "--cexcel", help="Complete Excel filename")
args = parser.parse_args()

with open(args.file) as inp:
    list = list(inp)
ID_list = []
for item in list:
    ID_list.append(item.replace("\n",""))

cexcel = Annotation.CompleteExcel(args.cexcel)
cexcel.getPeptidefa(ID_list, args.output)
