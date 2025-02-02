
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Fasta file")
parser.add_argument("-e", "--ematrix", help="Expression matrix")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()


miRNA_ids_file = "miRNA_ids.txt"  # Replace with the path to your miRNA IDs file
fasta_file = "miRNA_sequences.fa"  # Replace with the path to your FASTA file
output_file = "filtered_sequences.fa"  # Replace with the desired output file path

# Read miRNA IDs from the file
with open(args.ematrix, 'r') as ematrix:
    miRNA_id_list = [line.strip().split("\t")[0] for line in ematrix]

# Read the FASTA file and filter sequences based on miRNA IDs
with open(args.file, 'r') as fasta, open(args.output, 'wt') as output:
    current_id = None
    include_sequence = False

    for line in fasta:
        if line.startswith('>'):
            current_id = line[1:].strip()
            include_sequence = current_id in miRNA_id_list
        if include_sequence:
            output.write(line)
