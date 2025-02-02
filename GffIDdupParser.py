import argparse
import Annotation

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Gff file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


def main(input_file, output_file):
    ids = {}
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in sorted(infile, key=lambda x: (x.split('\t')[0], int(x.split('\t')[3]), float(x.split('\t')[5])),
                           reverse=True):
            modified_line = Annotation.process_line(line, ids)
            outfile.write(modified_line + '\n')


if __name__ == "__main__":
    main(args.file, args.output)