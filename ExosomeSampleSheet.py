#!/usr/bin/env/python3

import argparse
import Annotation
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="ematrix file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


def exosome_sample_sheet(file,output):
    em = Annotation.Ematrix(file)
    columns = em.header.replace("\n","").split("\t")[1::]
    files = []
    groups = []
    samples = []
    for column in columns:
        if column.startswith("Undetermined"):
            file = "Undetermined_S0"
            sample = "Undetermined"
            group = "Undetermined"
        else:
            file = column
            sample = column.split("_")[0]
            size = sample.split("-")[1]
            if sample.startswith("F"):
                group = "FTS_%s" %size
            elif sample.startswith("1"):
                group = "1EXP_%s" % size
            elif sample.startswith("2"):
                group = "2EXP_%s" % size
        files.append(file)
        groups.append(group)
        samples.append(sample)

    df = pd.DataFrame({"file":files, "group":groups, "sample":samples})
    df.to_csv(output, index=False, sep="\t")


exosome_sample_sheet(args.file, args.output)