from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Transcriptome fasta file")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

# Open the GTF file for writing
with open(args.output, "wt") as gtf:

    # Iterate through the FASTA sequences
    for record in SeqIO.parse(args.file, "fasta"):
        # Get sequence information
        seq_id = record.id
        seq_length = len(record.seq)

        # Define GTF fields
        seq_name = seq_id
        feature_type = "transcript"
        start = 1  # Start coordinate
        end = seq_length  # End coordinate
        score = "."
        strand = "+"
        frame = "."
        attributes = 'gene_id "{}";'.format(seq_id)

        # Write GTF entry
        gtf_line = "\t".join(
            [
                seq_name,
                "source",
                feature_type,
                str(start),
                str(end),
                score,
                strand,
                frame,
                attributes,
            ]
        )
        gtf.write(gtf_line + "\n")

