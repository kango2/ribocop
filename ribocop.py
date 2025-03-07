import pandas as pd
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt
import gzip
import tarfile
import json
import glob

def check_output(outputdir, sampleid):
    """Check if output has already been generated for a given sample ID."""
    checkpoint_file = os.path.join(outputdir, f"{sampleid}.ribocop")
    if os.path.exists(checkpoint_file + ".done"):
        print(f"Task already completed for {sampleid}. Exiting")
        sys.exit(0)
    log_file = os.path.join(outputdir, f"{sampleid}_log.json")
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            content = f.read()
            if 'Errors' in content.lower():
                print(f"{sampleid} previously produced error. Exiting")
                sys.exit(0)
    open(checkpoint_file + ".running", 'w').close()



def get_frequent_rdna(sampleid, log_file):

    """From a paf file, identify all the alignments with alignment length >90% of target length. identify the target 18s and 28s that occur most frequently in the input fasta (to be used as reference for second alingment)"""
    
    paf_file = f"{sampleid}.primary.paf"

    #redefine the column names in the paf file
    col_names = [
        "query_name", "query_length", "query_start", "query_end", 
        "strand", "target_name", "target_length", "target_start", "target_end", 
        "num_matches", "alignment_length", "mapping_quality", "tp", "cm", "s1", "s2", "dv", "rl"    ]

    # Read in PAF file and filter alignments with length > 90% of target length
    try:
        data = pd.read_csv(paf_file, sep='\t', names=col_names)
    except FileNotFoundError:
        print(f"Error: File {paf_file} not found!")
        sys.exit(1)   
 
    data["target_length"] = pd.to_numeric(data["target_length"], errors="coerce")
    data["alignment_length"] = pd.to_numeric(data["alignment_length"], errors="coerce")

    #Update log file with number of alignments
    update_log("rDNA_details", "Number of primary alignments", int(data.shape[0]), log_file)

    #Create directory to store figures in, if it does not already exist
    os.makedirs("figures", exist_ok=True)

    #Create an alignment length distribution plot
    plt.hist(data["alignment_length"], bins=30)
    plt.title(f"Distribution of alignment lengths for {sampleid}")
    plt.savefig(f"figures/{sampleid}.png", bbox_inches='tight')
    plt.clf()  # Clear the figure

    data = data[data["alignment_length"] > data["target_length"] * 0.85]

    #Update log file with number of alignments
    update_log("rDNA_details", "Number of filtered primary alignments", int(data.shape[0]), log_file)


    plt.hist(data["alignment_length"], bins=20)
    plt.title(f"Distribution of filtered alignment lengths for {sampleid}")
    plt.savefig(f"figures/{sampleid}_filtered.png", bbox_inches='tight')
    plt.clf()  # Clear the figure

    data["target_name"] = data["target_name"].astype(str)

    #Identify 18 and 28s target sequences with the most alignments
    filtered = data.loc[data["target_name"].str.contains("_18S", na=False), "target_name"].value_counts()
    eighteen = filtered.idxmax() if not filtered.empty else None

    total_18s_alignments = filtered.sum() if not filtered.empty else 0
    eighteen_count = filtered.max() if not filtered.empty else 0

    #Update log with details about 18S alignments
    update_log("rDNA_details", "Primary 18S", eighteen, log_file)
    update_log("rDNA_details", "Number of filtered primary 18S alignments", int(total_18s_alignments), log_file)
    update_log("rDNA_details", "Number of filtered primary alignments to chosen 18S", int(eighteen_count), log_file)

    if eighteen is None:
        print(f"No 18S found for {sampleid}. Exiting")
        update_log("Errors", "rDNA identification", f"No 18S found for {sampleid}.", log_file)
        exit()

    filtered = data.loc[data["target_name"].str.contains("_28S", na=False), "target_name"].value_counts()
    twoeight = filtered.idxmax() if not filtered.empty else None

    total_28s_alignments = filtered.sum() if not filtered.empty else 0
    twoeight_count = filtered.max() if not filtered.empty else 0

    update_log("rDNA_details", "Primary 28S", twoeight, log_file)
    update_log("rDNA_details", "Number of filtered primary 28S alignments", int(total_28s_alignments), log_file)
    update_log("rDNA_details", "Number of filtered primary alignments to chosen 28S", int(twoeight_count), log_file)

    if twoeight is None:
        print(f"No 28S found for {sampleid}. Exiting")
        update_log("Errors", "rDNA identification", f"No 28S found for {sampleid}.", log_file)
        exit()

    return eighteen, twoeight

def primary_alignment(rdnalibfa, PBS_NCPUS, sampleid, outputdir, inputfasta, log_file):
    paf_file = f"{sampleid}.primary.paf"
    #Align the input fasta to the rDNA library fasta. use get_frequent_rdna to identify the most commonly occuring 18s and 28s target sequences in the input fasta. extract these sequences from the rDNA library and write to a new fasta file
    
    command = f"minimap2 -t {PBS_NCPUS} --secondary=no -o {paf_file} {rdnalibfa} {inputfasta}"

    print(f"Running command: {command}")
    
    # Execute minimap2 command
    subprocess.run(command, check=True, shell=True)
    
    eighteen, twoeight = get_frequent_rdna(sampleid, log_file)
    print(eighteen, twoeight)

    with open(f'{outputdir}/{sampleid}.ssulsurdna.fa', 'w') as outfile:
        for record in SeqIO.parse(rdnalibfa, "fasta"):
        # Check if the record's ID matches one of the target IDs
            if record.id == eighteen or record.id == twoeight:
                SeqIO.write(record, outfile, 'fasta')


def refined_alignment(PBS_NCPUS, sampleid, inputfasta):
    #Align input fasta to our chosen refined rDNA sequences. 

    command = f"minimap2 -t {PBS_NCPUS} --secondary=no -o {sampleid}.refined.paf {sampleid}.ssulsurdna.fa {inputfasta}"

    print(f"Running command: {command}")

    subprocess.run(command, check=True, shell=True)


def read_paf_file(refined_paf):
    # Redefine the column names for the first 12 columns
    column_names = [
        "Query sequence name", "Query sequence length", "Query start", "Query end", "Relative strand",
        "Target sequence name", "Target sequence length", "Target start", "Target end", "Number of residue matches",
        "Alignment block length", "Mapping quality"
    ]

    # Reading the file with more control over parsing
    data = []
    with open(refined_paf, 'r') as file:
        for line in file:
            # Split the line by tab
            parts = line.strip().split('\t')
            data.append(parts)
    # Convert the list of lists to a DataFrame
    data_df = pd.DataFrame(data)
    # Assign column names to the first 12 columns
    data_df.columns = column_names + [f'Var_col_{i}' for i in range(1, len(data_df.columns) - 12 + 1)]
    # Convert relevant columns to numeric types
    numeric_columns = ["Query sequence length", "Query start", "Query end", 
                       "Target sequence length", "Target start", "Target end", "Number of residue matches",
                       "Alignment block length", "Mapping quality"]
    for col in numeric_columns:
        data_df[col] = pd.to_numeric(data_df[col])

    
    return data_df


def process_paf_alignments(data_df, sampleid, log_file):
    
    #Update log file
    update_log("rDNA_details", "Number of refined alignments", int(data_df.shape[0]), log_file)
    update_log("rDNA_details", "Median ratio of alignment block length to target sequence length for refined alignments", (data_df["Alignment block length"]/data_df["Target sequence length"]).median(), log_file)
    update_log("rDNA_details", "Median ratio of residue matches to target sequence length for refined alignments", (data_df["Number of residue matches"]/data_df["Target sequence length"]).median(), log_file)

    # Filter out short alignments, retaining alignments that cover at least 95% of the target sequence
    data_df = data_df[data_df['Alignment block length'].astype(float) / data_df['Target sequence length'].astype(float) >= 0.85]
    data_df = data_df[data_df['Number of residue matches'].astype(float) / data_df['Target sequence length'].astype(float) >= 0.50]
    #Previous results (in chordata directory) have 80% ratio requirement for number of residue matches 

    #Update log file
    update_log("rDNA_details", "Number of filtered refined alignments", int(data_df.shape[0]), log_file)
    update_log("rDNA_details", "Median ratio of alignment block length to target sequence length for filtered refined alignments", (data_df["Alignment block length"]/data_df["Target sequence length"]).median(), log_file)
    update_log("rDNA_details", "Median ratio of residue matches to target sequence length for filtered refined alignments", (data_df["Number of residue matches"]/data_df["Target sequence length"]).median(), log_file)


   # data_df = data_df.groupby("Target sequence name")
   # data_df = data_df[numpy.percentile(data_df.'Number of residue matches', 95)]

    #percentile_20 = data_df.groupby("Target sequence name")["Number of residue matches"].transform(lambda x: np.percentile(x, 20))

    # Filter rows where 'Number of residue matches' is above the 90th percentile
    #data_df = data_df[data_df["Number of residue matches"] >= percentile_20]

    # Ensure 18S is followed by 28S or 28S is followed by 18S for each query sequence name
    # Current logic captures 28S-18S units. we need to change it to 18S-28S units.
    morphs = set()
    for query_name, group in data_df.groupby("Query sequence name"):
        group = group.sort_values(by=["Query start"]).reset_index(drop=True)
        valid_group = set()
        i = 0
        while i < len(group):
            # start with a 18S in the plus strand, and 100bp available before this 18S on the contig
            if '_18S' in group.at[i, 'Target sequence name'] and group.at[i, 'Relative strand'] == "+" and group.at[i, 'Query start'] - 100 > 0:
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '_28S' in group.at[i + 1, 'Target sequence name'] and group.at[i + 1, 'Relative strand'] == "+":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '_18S' in group.at[i + 2, 'Target sequence name'] and group.at[i + 2, 'Relative strand'] == "+":
                        morph = group.at[i, 'Query sequence name'], group.at[i, 'Query start'] - 100, group.at[i + 2, 'Query start'] - 100, '+'
                        morphs.add(morph)
            # start with a 18S in the minus strand
            elif '_18S' in group.at[i, 'Target sequence name'] and group.at[i, 'Relative strand'] == "-":
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '_28S' in group.at[i + 1, 'Target sequence name'] and group.at[i, 'Relative strand'] == "-":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '_18S' in group.at[i + 2, 'Target sequence name'] and group.at[i + 2, 'Relative strand'] == "-" and group.at[i + 2, 'Query end'] + 100 <= group.at[i, 'Query sequence length']:
                        morph = group.at[i, 'Query sequence name'], group.at[i, 'Query end'] + 100, group.at[i + 2, 'Query end'] + 100, '-'
                        morphs.add(morph)                    

            i += 1
    # Convert morphs set to DataFrame
    if len(morphs)==0:
        print(f"No morphs pass filtering for {sampleid}.Exiting")
        update_log("Errors", "Morph identification", f"No morphs pass filtering for {sampleid}.", log_file)
        exit()
    else:
        morphs_df = pd.DataFrame(list(morphs), columns=['Query sequence name', 'Start', 'End', 'Strand'])
        update_log("rDNA_details", "Number of morphs", len(morphs), log_file)
    return morphs_df
    # return morphs

def extract_sequences(fasta_file, morphs_df, output_dir, sampleid):
    if fasta_file.endswith(".gz"):
        with gzip.open(fasta_file, "rt") as handle:
            sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        with open(fasta_file, "rt") as handle:
            sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    extracted_sequences = []
    valid_rows = []

    for _, row in morphs_df.iterrows():
        seq_id = row['Query sequence name']
        start = int(row['Start'])
        end = int(row['End'])
        strand = row['Strand']

        if seq_id in sequences:
            sequence = sequences[seq_id].seq[start:end]
            if strand == '-':
                sequence = sequence.reverse_complement()
            sequence_str = str(sequence)
            if 'N' not in sequence_str:
                extracted_sequences.append((seq_id + ':' + str(start) + '-' + str(end) + ':' + strand, sequence_str))
                valid_rows.append(row)

    with open(os.path.join(output_dir, sampleid + ".rDNA.morphs.fasta"), 'w') as output_handle:
        for seq_id, sequence in extracted_sequences:
            output_handle.write(f">{seq_id}\n{sequence}\n")

    valid_morphs_df = pd.DataFrame(valid_rows)
    output_file = os.path.join(output_dir, sampleid + ".rDNA.morphs.tsv")
    valid_morphs_df.to_csv(output_file, sep='\t', index=False, header=False)

def rna_builder(sampleid, inputfasta, outputdir, log_file):
    refined_paf = f"{sampleid}.refined.paf"

    data_df = read_paf_file(refined_paf)
    morphs_df = process_paf_alignments(data_df, sampleid, log_file)
    extract_sequences(inputfasta, morphs_df, outputdir, sampleid)

def get_median_morph(sampleid, PBS_NCPUS, inputfasta, log_file):
    morph_fasta = f"{sampleid}.rDNA.morphs.fasta"

    fai_file = f"{morph_fasta}.fai"

    subprocess.run(f"samtools faidx {morph_fasta}", shell=True, check=True)

    morphs = pd.read_csv(fai_file, sep='\t', header=None, names=['sequence_name', 'length', 'offset', 'linebases', 'linewidth', 'qualoffset'])

    # Sort by length
    morphs_sorted = morphs.sort_values(by='length')

    # Find the median length sequence
    median_index = len(morphs_sorted) // 2
    refrdnamorph = morphs_sorted.iloc[median_index]['sequence_name']
    refmorphlength = morphs_sorted.iloc[median_index]['length']
    min_length = morphs_sorted.iloc[0]['length']
    max_length = morphs_sorted.iloc[-1]['length']


    #Update log file
    update_log("rDNA_details", "Unit length", int(refmorphlength), log_file)
    update_log("rDNA_details", "Minimum unit length", int(min_length), log_file)
    update_log("rDNA_details", "Maximum unit length", int(max_length), log_file)

    command = f"samtools faidx {morph_fasta} {refrdnamorph} > {sampleid}.rDNA.refmorph.fasta"

    print(f"Running command: {command}")

    subprocess.run(command, check=True, shell=True)


    command=f"minimap2 -t {PBS_NCPUS} --secondary=no -o {sampleid}.asm2refmorph.paf {sampleid}.rDNA.refmorph.fasta {inputfasta}"
    print(f"Running command: {command}")

    subprocess.run(command, shell=True, check=True)

    

def process_paf(sampleid, outputdir, log_file):

    column_names = [
        "Query sequence name", "Query sequence length", "Query start", "Query end", "Relative strand",
        "Target sequence name", "Target sequence length", "Target start", "Target end", "Number of residue matches",
        "Alignment block length", "Mapping quality"
    ]

    #Read in paf file as pandas df and set output file
    data = pd.read_csv(f"{sampleid}.asm2refmorph.paf", sep='\t', usecols=list(range(0,12)), names=column_names)
    outputfile = f"{sampleid}.rdnaregions.bed"

    #Update log file
    update_log("rDNA_details", "Number of final alignments", int(data.shape[0]), log_file)
 

    #Set variables
    alnratio = 0.8
    minalnresidue = 500
    edge = 100
    gapwidth = 100

    filtered_data = data[
        (data["Number of residue matches"].astype(float) > minalnresidue) & 
        (data["Number of residue matches"].astype(float) / data["Alignment block length"].astype(float) > alnratio)
    ]

    #plt.hist(filtered_data.loc[filtered_data["Target sequence name"].str.contains("_18S", na=False), "Alignment block length"], bins=20)
    #plt.title(f"Filtered and refined morph lengths for {sampleid} 18S")
   # plt.savefig(f"figures/{sampleid}_final_18S.png", bbox_inches='tight')
   # plt.clf()  # Clear the figure

    #plt.hist(filtered_data.loc[filtered_data["Target sequence name"].str.contains("_28S", na=False), "Alignment block length"], bins=20)
    ##plt.title(f"Filtered and refined morph lengths for {sampleid} 28S")
    #plt.savefig(f"figures/{sampleid}_final_28S.png", bbox_inches='tight')
   # plt.clf()  # Clear the figure

    minregionlength = filtered_data["Number of residue matches"].median() / 2

    filtered_data.loc[:, "Query start"] = np.where(
        filtered_data["Query start"] < edge, 
        0, 
        filtered_data["Query start"]
    )

    filtered_data.loc[:, "Query end"] = np.where(
        (filtered_data["Query sequence length"] - filtered_data["Query end"]) < edge, 
        filtered_data["Query sequence length"], 
        filtered_data["Query end"]
    )

    filtered_data = filtered_data.sort_values(["Query sequence name", "Query start"]).reset_index(drop=True)

    filtered_data["group"] = (filtered_data["Query start"] > filtered_data.groupby("Query sequence name")["Query end"].cummax().shift(1).fillna(-float("inf")) + gapwidth).cumsum()

    result = filtered_data.groupby(["Query sequence name", "group"]).agg({"Query start": "min", "Query end": "max"}).reset_index()

    result = result[["Query sequence name", "Query start", "Query end"]]

    result = result[(result["Query end"] - result["Query start"]) > minregionlength]

    total_length = (result['Query end'] - result['Query start']).sum()

    contigs = len(pd.unique(result["Query sequence name"]))

    #Update log file
    update_log("rDNA_details", "Total length", int(total_length), log_file)
    update_log("rDNA_details", "Number of contigs", int(contigs), log_file)
    update_log("rDNA_details", "Number of final filtered alignments", int(result.shape[0]), log_file)

    result.to_csv(outputfile, sep='\t', index=False, header=False)



        
def update_log(category, key, value, log_file):
    with open(log_file, "r") as f:
            try:
                data = json.load(f)  # Load existing JSON data
            except json.JSONDecodeError:
                data = {}  # Handle corrupted/empty file
    if category not in data:
        data[category]={}

    data[category][key] = value

    with open(log_file, "w") as f:
        json.dump(data, f, indent=2)


def main():
    """Main function."""   
    parser = argparse.ArgumentParser(description="Identify and characterise rDNA units from FASTA input.")
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for results')
    parser.add_argument('-s', '--sampleid', required=True, help='Sample ID for output files')
    parser.add_argument('-l', '--rdnalib', required=True, help='rDNA library')
    parser.add_argument('-t', '--ncpus', required=True, help='Threads')
    parser.add_argument('-x', '--txid', required=True, help='txID')
    parser.add_argument('-p', '--species', required=True, help='Species name')
    parser.add_argument('-i', '--input_dir', required=False, default=".", help="Input directory (default: current directory)")



    args = parser.parse_args()

    
#Create output directory if it does not exist.
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"Directory '{args.output_dir}' created successfully.")
    except OSError as e:
        print(f"Error creating directory '{args.output_dir}': {e}")
        sys.exit(1)

#Check if output already exists
    check_output(args.output_dir, args.sampleid)

    #Set wd to output directory
    try:
        os.chdir(args.output_dir)
        print("WD changed successfully")
    except Exception as e:
        print(f"Error changing directory: {e}")
        sys.exit(1)

    log_file = os.path.join(args.output_dir, f"{args.sampleid}_log.json")

    sampleinfo = {
        "Sample_details":{
        "Sample ID": args.sampleid,
        "TaxID":args.txid,
        "Species":args.species
        }
    }

    with open(log_file, "w") as f:
        json.dump(sampleinfo, f, indent=2)

    inputfasta = (
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fna") + 
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fna.gz")
    )[0]
    
    primary_alignment(args.rdnalib, args.ncpus, args.sampleid, args.output_dir, inputfasta, log_file)
    refined_alignment(args.ncpus, args.sampleid, inputfasta)

    rna_builder(args.sampleid, inputfasta, args.output_dir, log_file)

    get_median_morph(args.sampleid, args.ncpus, inputfasta, log_file)

    process_paf(args.sampleid, args.output_dir, log_file)

    #summary_file = os.path.join(args.output_dir, "summary.json")

 

    checkpoint_file = os.path.join(args.output_dir, f"{args.sampleid}.ribocop")
    
    open(checkpoint_file + ".done", 'w').close()

    os.remove(f"{checkpoint_file}.running")

    files = [f for f in os.listdir() if f.startswith(f"{args.sampleid}")]
    with tarfile.open(name=f"{args.sampleid}.tar.gz", mode="w:gz") as tar:
        for file in files:
            tar.add(file)

    print("Task completed. Exiting")
    
    

if __name__ == "__main__":
    main()     
