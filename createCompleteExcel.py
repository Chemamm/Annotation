#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--excel", help="excel best alignment")
parser.add_argument("-o", "--output", help="output filename")
parser.add_argument("-f", "--fasta", help="fasta transcriptome filename")
parser.add_argument("-c", "--cds", help="cds fasta filename")
parser.add_argument("-p", "--pep", help="cds pep filename")
parser.add_argument("-m", "--matrix", help="fpkm expresion matrix")
parser.add_argument("-u", "--uniprot", help="uniprot annotation file")
parser.add_argument("-i", "--interpro", help="interpro annotation file")
parser.add_argument("-b", "--blast2go", help="blast2go annotation file")
parser.add_argument("-s", "--signalp", help="signalp annotation file")
parser.add_argument("-ph", "--phobius", help="phobius annotation file")
parser.add_argument("-t", "--tmhmm", help="tmhmm annotation file")
parser.add_argument("-tm", "--tmhmm_mature", help="tmhmm mature annotation file")
parser.add_argument("-go", "--go_obo", help="GO anottation obo file")
parser.add_argument("-de", "--list_of_DEfiles", help="list of DE files. Two columns: filename and code")
args = parser.parse_args()

Annotation.createCompleteExcel(args.excel, args.fasta, args.cds, args.pep, args.matrix, args.uniprot, args.blast2go, args.interpro, args.signalp, args.phobius, args.tmhmm_mature, args.tmhmm, args.output, args.go_obo, args.list_of_DEfiles)
