import pandas as pd
import gzip
import os
import argparse
import pathlib
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="Expression Matrix")
parser.add_argument("-fq", "--fastq_dir", help="Fastq folder")
parser.add_argument("-o", "--output", help="output filename")
args = parser.parse_args()

# Dictionary to map old sample names to new sample names
sample_name_mapping = {
    'pre_2_SG_S235_L004_Aligned.sortedByCoord.out.bam': 'SG_unfed_pre2',
    'pre_4_SG_S267_L004_Aligned.sortedByCoord.out.bam': 'SG_unfed_pre4',
    'pre_6_SG_S236_L004_Aligned.sortedByCoord.out.bam': 'SG_unfed_pre6',
    'pre_8_SG_S274_L004_Aligned.sortedByCoord.out.bam': 'SG_unfed_pre8',
    'pre_10_SG_S237_L004_Aligned.sortedByCoord.out.bam': 'SG_unfed_pre10',
    '1_12h-SG_S268_L004_Aligned.sortedByCoord.out.bam': 'SG_1_12h_1',
    '105-12h-SG_S1_L001_Aligned.sortedByCoord.out.bam': 'SG_1_12h_105',
    '5_12h_SG_S224_L004_Aligned.sortedByCoord.out.bam': 'SG_1_12h_5',
    '101_12h_SG_S225_L004_Aligned.sortedByCoord.out.bam': 'SG_1_12h_101',
    '103_12h_SG_S275_L004_Aligned.sortedByCoord.out.bam': 'SG_1_12h_103',
    '7_24h_SG_S270_L004_Aligned.sortedByCoord.out.bam': 'SG_1_24h_7',
    '9_24h_SG_S226_L004_Aligned.sortedByCoord.out.bam': 'SG_1_24h_9',
    '109_24h_SG_S312_L004_Aligned.sortedByCoord.out.bam': 'SG_1_24h_109',
    '106_24h_SG_S227_L004_Aligned.sortedByCoord.out.bam': 'SG_1_24h_106',
    '8-24h-SG_S269_L004_Aligned.sortedByCoord.out.bam': 'SG_1_24h_8',
    '13_48h_SG_S255_L004_Aligned.sortedByCoord.out.bam': 'SG_1_48h_13',
    '15_48h_SG_S256_L004_Aligned.sortedByCoord.out.bam': 'SG_1_48h_15',
    '111_48h_SG_S228_L004_Aligned.sortedByCoord.out.bam': 'SG_1_48h_111',
    '113_48h_SG_S229_L004_Aligned.sortedByCoord.out.bam': 'SG_1_48h_113',
    '114_48h_SG_S271_L004_Aligned.sortedByCoord.out.bam': 'SG_1_48h_114',
    '18_72h_SG_S257_L004_Aligned.sortedByCoord.out.bam': 'SG_1_72h_18',
    '19_72h_SG_S258_L004_Aligned.sortedByCoord.out.bam': 'SG_1_72h_19',
    '20_72h_SG_S259_L004_Aligned.sortedByCoord.out.bam': 'SG_1_72h_20',
    '116_72h_SG_S307_L004_Aligned.sortedByCoord.out.bam': 'SG_1_72h_116',
    '118_72h_SG_S230_L004_Aligned.sortedByCoord.out.bam': 'SG_1_72h_118',
    '21_96h_SG_S239_L004_Aligned.sortedByCoord.out.bam': 'SG_1_96h_21',
    '22_96h_SG_S240_L004_Aligned.sortedByCoord.out.bam': 'SG_1_96h_22',
    '24_96h_SG_S260_L004_Aligned.sortedByCoord.out.bam': 'SG_1_96h_24',
    '121_96h_SG_S231_L004_Aligned.sortedByCoord.out.bam': 'SG_1_96h_121',
    '122_96h_SG_S261_L004_Aligned.sortedByCoord.out.bam': 'SG_1_96h_122',
    '52_12h_SG_S280_L004_Aligned.sortedByCoord.out.bam': 'SG_2_12h_52',
    '53_12h_SG_S232_L004_Aligned.sortedByCoord.out.bam': 'SG_2_12h_53',
    '152_12h_SG_S233_L004_Aligned.sortedByCoord.out.bam': 'SG_2_12h_152',
    '153_12h_SG_S234_L004_Aligned.sortedByCoord.out.bam': 'SG_2_12h_153',
    '155_12h_SG_S272_L004_Aligned.sortedByCoord.out.bam': 'SG_2_12h_155',
    '56_1d_SG_S262_L004_Aligned.sortedByCoord.out.bam': 'SG_2_24h_56',
    '59_1d_SG_S263_L004_Aligned.sortedByCoord.out.bam': 'SG_2_24h_59',
    '156_1d_SG_S241_L004_Aligned.sortedByCoord.out.bam': 'SG_2_24h_156',
    '157_1d_SG_S273_L004_Aligned.sortedByCoord.out.bam': 'SG_2_24h_157',
    '159_1d_SG_S242_L004_Aligned.sortedByCoord.out.bam': 'SG_2_24h_159',
    '62_2d_SG_S243_L004_Aligned.sortedByCoord.out.bam': 'SG_2_48h_62',
    '63_2d_SG_S244_L004_Aligned.sortedByCoord.out.bam': 'SG_2_48h_63',
    '65_2d_SG_S245_L004_Aligned.sortedByCoord.out.bam': 'SG_2_48h_65',
    '162_2d_SG_S246_L004_Aligned.sortedByCoord.out.bam': 'SG_2_48h_162',
    '163_2d_SG_S264_L004_Aligned.sortedByCoord.out.bam': 'SG_2_48h_163',
    '67_3d_SG_S247_L004_Aligned.sortedByCoord.out.bam': 'SG_2_72h_67',
    '68_3d_SG_S265_L004_Aligned.sortedByCoord.out.bam': 'SG_2_72h_68',
    '69_3d_SG_S248_L004_Aligned.sortedByCoord.out.bam': 'SG_2_72h_69',
    '166_3d_SG_S249_L004_Aligned.sortedByCoord.out.bam': 'SG_2_72h_166',
    '167_3d_SG_S250_L004_Aligned.sortedByCoord.out.bam': 'SG_2_72h_167',
    '71_4d_SG_S266_L004_Aligned.sortedByCoord.out.bam': 'SG_2_96h_71',
    '74_4d_SG_S251_L004_Aligned.sortedByCoord.out.bam': 'SG_2_96h_74',
    '75_4d_SG_S252_L004_Aligned.sortedByCoord.out.bam': 'SG_2_96h_75',
    '171_4d_SG_S253_L004_Aligned.sortedByCoord.out.bam': 'SG_2_96h_171',
    '172_4d_SG_S254_L004_Aligned.sortedByCoord.out.bam': 'SG_2_96h_172',
    'pre4_MG_S24_L001_Aligned.sortedByCoord.out.bam': 'MG_unfed_pre4',
    'pre10-MG_S2_L001_Aligned.sortedByCoord.out.bam': 'MG_unfed_pre10',
    'pre_8_MG_S279_L004_Aligned.sortedByCoord.out.bam': 'MG_unfed_pre8',
    '4_12h_MG_S310_L004_Aligned.sortedByCoord.out.bam': 'MG_1_12h_4',
    '5_12h_MG_S311_L004_Aligned.sortedByCoord.out.bam': 'MG_1_12h_5',
    '104_12h_MG_S238_L004_Aligned.sortedByCoord.out.bam': 'MG_1_12h_104',
    '9_24h_MG_S297_L004_Aligned.sortedByCoord.out.bam': 'MG_1_24h_9',
    '10_24h_MG_S298_L004_Aligned.sortedByCoord.out.bam': 'MG_1_24h_10',
    '107_24h_MG_S299_L004_Aligned.sortedByCoord.out.bam': 'MG_1_24h_107',
    '13_48h_MG_S281_L004_Aligned.sortedByCoord.out.bam': 'MG_1_48h_13',
    '15_48h_MG_S282_L004_Aligned.sortedByCoord.out.bam': 'MG_1_48h_15',
    '114_48h_MG_S283_L004_Aligned.sortedByCoord.out.bam': 'MG_1_48h_114',
    '20_72h_MG_S300_L004_Aligned.sortedByCoord.out.bam': 'MG_1_72h_20',
    '116_72h_MG_S284_L004_Aligned.sortedByCoord.out.bam': 'MG_1_72h_116',
    '118_72h_MG_S285_L004_Aligned.sortedByCoord.out.bam': 'MG_1_72h_118',
    '22_96h_MG_S286_L004_Aligned.sortedByCoord.out.bam': 'MG_1_96h_22',
    '24_96h_MG_S287_L004_Aligned.sortedByCoord.out.bam': 'MG_1_96h_24',
    '122_96h_MG_S288_L004_Aligned.sortedByCoord.out.bam': 'MG_1_96h_122',
    '155_12h_MG_S308_L004_Aligned.sortedByCoord.out.bam': 'MG_2_12h_155',
    '54_12_MG_S309_L004_Aligned.sortedByCoord.out.bam': 'MG_2_12h_54',
    '152_12h_MG_S276_L004_Aligned.sortedByCoord.out.bam': 'MG_2_12h_152',
    '56_1d_MG_S277_L004_Aligned.sortedByCoord.out.bam': 'MG_2_24h_56',
    '59_1d_MG_S301_L004_Aligned.sortedByCoord.out.bam': 'MG_2_24h_59',
    '157_1d_MG_S278_L004_Aligned.sortedByCoord.out.bam': 'MG_2_24h_157',
    '62_2d_MG_S289_L004_Aligned.sortedByCoord.out.bam': 'MG_2_48h_62',
    '63_2d_MG_S290_L004_Aligned.sortedByCoord.out.bam': 'MG_2_48h_63',
    '163_2d_MG_S302_L004_Aligned.sortedByCoord.out.bam': 'MG_2_48h_163',
    '67_3d_MG_S291_L004_Aligned.sortedByCoord.out.bam': 'MG_2_72h_67',
    '68_3d_MG_S292_L004_Aligned.sortedByCoord.out.bam': 'MG_2_72h_68',
    '167_3d_MG_S293_L004_Aligned.sortedByCoord.out.bam': 'MG_2_72h_167',
    '71_4d_MG_S294_L004_Aligned.sortedByCoord.out.bam': 'MG_2_96h_71',
    '74_4d-MG_S295_L004_Aligned.sortedByCoord.out.bam': 'MG_2_96h_74',
    '171_4d_MG_S296_L004_Aligned.sortedByCoord.out.bam': 'MG_2_96h_171'
}


# Read the expression matrix
expression_df = pd.read_csv(args.file, sep='\t', index_col=0, comment='#')


# Define a function to count sequences in gzipped FASTQ files
def count_sequences_in_fastq(fastq_file):
    with gzip.open(fastq_file, 'rt') as f:
        return sum(1 for line_number, line in enumerate(f) if line_number % 4 == 1)


# File to store total reads counts
total_reads_file = '/shared/raw_data/Tick_Transcriptome/trimmed/read_counts.txt'

# Check if the total reads file exists
if pathlib.Path(total_reads_file).is_file():
    # If the file exists, read total reads counts from it
    with open(total_reads_file, 'r') as trf:
        total_reads_dict = {line.split()[0]: int(line.split()[1]) for line in trf}

else:
    # If the file doesn't exist, calculate total reads and save them to the file
    total_reads_dict = {}

    # Iterate through the columns of the expression matrix (excluding non-sample columns)
    for column in expression_df.columns[5:]:
        # Generate the corresponding FASTQ file name
        fastq_file = os.path.join(args.fastq_dir, column.replace('_Aligned.sortedByCoord.out.bam', '_R1_001_val_1.fq.gz'))

        # Count the sequences in the FASTQ file
        total_reads = count_sequences_in_fastq(fastq_file)

        # Store the total reads count in the dictionary
        total_reads_dict[column] = total_reads

    # Save the total reads dictionary to the file
    with open(total_reads_file, 'wt') as trf:
        for sample, total_reads in total_reads_dict.items():
            trf.write(f"{sample}\t{total_reads}\n")

order =  list(sample_name_mapping.values())

# Normalize the expression matrix to FPKMs
for column in expression_df.columns[5:]:
    sample_name = column
    if sample_name in sample_name_mapping:
        new_sample_name = sample_name_mapping[sample_name]
        expression_df.rename(columns={sample_name: new_sample_name}, inplace=True)
        total_reads_dict[new_sample_name] = total_reads_dict.pop(sample_name)

    total_reads = total_reads_dict[new_sample_name]
    expression_df[new_sample_name] = (expression_df[new_sample_name] / expression_df['Length']) / (total_reads / 1e9)

expression_df = expression_df[order]

# Save the FPKM matrix to a new file
expression_df.to_csv(args.output)

#Logarithmic values
log_df = np.log(expression_df)
log_df.to_csv(args.output.replace(".txt","_log.txt"))
