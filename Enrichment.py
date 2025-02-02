#!/usr/bin/env/python3

import argparse
import Annotation
import os


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--list_of_files", help="list of DE files and its code")
parser.add_argument("-e", "--excel", help="Complete Excel")
parser.add_argument("-o", "--output", help="background filescoutput label")
parser.add_argument("-p", "--pwd", help="path to working directory")
args = parser.parse_args()

os.chdir(args.pwd)
ex = Annotation.CompleteExcel(args.excel)
ex.getBackgrounds(args.output)
filelist = Annotation.DEgetlist(args.list_of_files, ex.transcriptIDs)
Annotation.enrichment(filelist, args.output, args.pwd)
