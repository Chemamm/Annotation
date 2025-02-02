#!/usr/bin/env/python3

import argparse
import Annotation


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Complete Excel")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

ce = Annotation.CompleteExcel(args.file)

with open(args.output, "wt") as out:
    out.write("\n".join(ce.secreted))
