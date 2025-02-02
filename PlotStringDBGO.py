#!/usr/bin/env/python3

import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-bp", "--biological_process", default=False, help="biological process GO term file")
parser.add_argument("-mf", "--molecular_function", default=False, help="molecular function GO term file")
parser.add_argument("-cc", "--celular_component", default=False, help="celular component GO term file")
parser.add_argument("-cl", "--custom_list", default=False, help="custom list GO term file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


Annotation.stringDB_enrichment_figure(args.biological_process,args.molecular_function,args.celular_component,
                                      args.custom_list, args.output)

