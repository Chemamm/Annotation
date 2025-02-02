import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
import os

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Select best protein alignments and generate ClustalW alignments.")
    parser.add_argument("-b", "--blast_result", required=True, help="BLAST result file in outfmt 6 format.")
    parser.add_argument("-q", "--query_fasta", required=True, help="Query protein sequences in FASTA format.")
    parser.add_argument("-s", "--subject_fasta", required=True, help="Subject protein sequences in FASTA format.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for ClustalW alignments.")
    parser.add_argument("-f", "--final_summary", required=True, help="Output TSV file for final summary of best alignments.")
    return parser.parse_args()

def select_best_hits(blast_result_file):
    """Parse BLAST results and select the best alignment for each query sequence."""
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qlen",
               "qstart", "qend", "slen", "sstart", "send", "evalue", "bitscore"]
    blast_df = pd.read_csv(blast_result_file, sep="\t", names=columns)

    # Sort by query ID and bit score, then keep only the best hit per query
    best_hits = blast_df.sort_values(by=["qseqid", "bitscore"], ascending=[True, False]).drop_duplicates("qseqid")
    return best_hits

def extract_sequences(fasta_file):
    """Load sequences from a FASTA file into a dictionary."""
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

def align_sequences(query_seq, subject_seq, output_file):
    """Align two sequences using Clustal Omega and save the output in ClustalW format."""
    temp_fasta = output_file.replace(".clustalw", ".temp.fasta")

    # Write query and subject sequences to a temporary FASTA file
    with open(temp_fasta, "w") as temp_fasta_handle:
        SeqIO.write(query_seq, temp_fasta_handle, "fasta")
        SeqIO.write(subject_seq, temp_fasta_handle, "fasta")

    # Generate alignment in ClustalW format, allowing overwriting with --force
    clustalomega_cline = ClustalOmegaCommandline(
        infile=temp_fasta,
        outfile=output_file,
        verbose=True,
        auto=True,
        force=True,
        outfmt="clu"  # Add this flag to force overwrite
    )
    clustalomega_cline()

    # Clean up the temporary FASTA file
    os.remove(temp_fasta)

def calculate_additional_metrics(row):
    """Calculate query and subject coverage, and amino acids lacking at the end of the subject."""
    qseq_cover = (row['length'] / row['qlen']) * 100  # Query coverage percentage
    sseq_cover = (row['length'] / row['slen']) * 100  # Subject coverage percentage
    aas_lacking = row['slen'] - max(row['send'], row['sstart'])  # Amino acids not aligning at the end
    return qseq_cover, sseq_cover, aas_lacking

def main():
    args = parse_arguments()

    # Step 1: Parse BLAST results
    print("Parsing BLAST results...")
    best_hits = select_best_hits(args.blast_result)

    # Step 2: Load query and subject sequences
    print("Loading sequences...")
    query_sequences = extract_sequences(args.query_fasta)
    subject_sequences = extract_sequences(args.subject_fasta)

    # Step 3: Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Step 4: Generate alignments for the best hits
    print("Generating alignments...")
    summary_data = []  # Store data for the summary TSV file

    for _, row in best_hits.iterrows():
        query_id = row["qseqid"]
        subject_id = row["sseqid"]

        if query_id in query_sequences and subject_id in subject_sequences:
            query_seq = query_sequences[query_id]
            subject_seq = subject_sequences[subject_id]

            output_file = os.path.join(args.output_dir, f"{query_id}_alignment.clustalw")
            align_sequences(query_seq, subject_seq, output_file)
            print(f"Alignment saved: {output_file}")

            # Calculate additional metrics
            qseq_cover, sseq_cover, aas_lacking = calculate_additional_metrics(row)
            row_data = row.to_dict()
            row_data.update({
                "qseq_cover": round(qseq_cover, 2),
                "sseq_cover": round(sseq_cover, 2),
                "aas_lacking": int(aas_lacking)
            })
            summary_data.append(row_data)
        else:
            print(f"Warning: Sequences for {query_id} or {subject_id} not found!")

    # Step 5: Write final summary TSV file
    print("Writing final summary file...")
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(args.final_summary, sep="\t", index=False)
    print(f"Summary file saved: {args.final_summary}")

    print("All alignments completed successfully.")

if __name__ == "__main__":
    main()
