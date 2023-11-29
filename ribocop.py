import sys
import argparse
import pysam
import pandas as pd
import matplotlib.pyplot as plt

def read_fasta_file(fasta_file):
    sequence_lengths = {}
    with pysam.FastaFile(fasta_file) as fasta:
        for seq in fasta.references:
            sequence_lengths[seq] = fasta.get_reference_length(seq)
    return sequence_lengths

def read_alignment_file(filename, ref_lengths, unique_sequences):
    alnInfo = []  # List to store data for DataFrame

    # Open the SAM/BAM/CRAM file
    with pysam.AlignmentFile(filename, "r") as file:
        for alignment in file:
            # Check if the alignment is mapped
            if not alignment.is_unmapped:
                # Extract bases aligned to the reference
                # aligned_bases = alignment.query_sequence
                # Extract read length
                # read_length = alignment.query_length

                # Calculate the length of the aligned bases
                aligned_bases_length = alignment.query_alignment_length
                if aligned_bases_length < 500:
                    continue
                # Extract alignment score from AS tag
                alignment_score = alignment.get_tag('AS') if alignment.has_tag('AS') else None
                if alignment_score < 300:
                    continue
                # Get reference length from the dictionary
                ref_length = ref_lengths.get(alignment.reference_name, None)
                # normalised score
                normalised_aln_score = alignment_score / ref_length
                # Split reference name
                refid, refType = alignment.reference_name.split('_') if '_' in alignment.reference_name else (alignment.reference_name, None)

                # Insert query sequence into the dictionary
                if alignment.query_name not in unique_sequences:
                    unique_sequences[alignment.query_name] = alignment.query_sequence

                # Add data to the list
                alnInfo.append({
                    "query_name": alignment.query_name,
                    "reference_name": alignment.reference_name,
                    "ref_length": ref_length,
                    "aligned_bases": aligned_bases_length,
                    "normalised_aln_score": normalised_aln_score,
                    "refid": refid,
                    "refType": refType
                })
                # # Print the results for each alignment
                # print(f"Aligned Bases: {aligned_bases_length}")
                # print(f"Alignment Score: {alignment_score}")
                # print(f"Score per Base: {score_per_base}")
                # print(f"rlen:", ref_lengths.get(alignment.reference_name, None))
                # #print(f"Read Length: {read_length}")
                
                # print("-" * 40)
    return pd.DataFrame(alnInfo)

def process_bam_files(bam_files, output_prefix, ref_lengths):
    unique_sequences = {}
    combined_alnInfo = pd.DataFrame()

    for bam_file in bam_files:
        alnInfo = read_alignment_file(bam_file, ref_lengths, unique_sequences)
        alnInfo['filename'] = bam_file
        combined_alnInfo = pd.concat([combined_alnInfo, alnInfo])

    combined_alnInfo.to_csv(f'{output_prefix}_alignment_data.csv', index=False)

    with open(f'{output_prefix}_unique_sequences.fasta', 'w') as file:
        for query_name, sequence in unique_sequences.items():
            file.write(f">{query_name}\n{sequence}\n")

    for refType, group in combined_alnInfo.groupby('refType'):
        plt.figure()
        group['normalised_aln_score'].hist(bins=100)
        plt.title(f"Histogram of normalised_aln_score for refType: {refType}")
        plt.xlabel('normalised_aln_score')
        plt.ylabel('Frequency')
        plt.savefig(f"{output_prefix}_{refType}_histogram.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM files and generate alignment data.")
    parser.add_argument('-p', '--pacbiobam', nargs='+', help='BAM files generated using PacBio data')
    parser.add_argument('-o', '--ontbam', nargs='+', help='BAM files generated using ONT data')
    parser.add_argument('-r', '--reference', required=True, help='Reference FASTA file')
    parser.add_argument('-d', '--output_dir', default='.', help='Output directory for the files')
    args = parser.parse_args()

    ref_lengths = read_fasta_file(args.reference)

    if args.pacbiobam:
        process_bam_files(args.pacbiobam, f'{args.output_dir}/pacbio', ref_lengths)

    if args.ontbam:
        process_bam_files(args.ontbam, f'{args.output_dir}/ont', ref_lengths)

