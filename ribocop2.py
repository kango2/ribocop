import pandas as pd
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import sys
import numpy as np
import gzip
import tarfile
import json
import glob
import re


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




def run_hmm(PBS_NCPUS, sampleid, outputdir, inputfasta, log_file, hmmdb):
    #Align the input fasta to HMM library using nhmmer.
    PBS_JOBFS = os.environ.get('PBS_JOBFS')

    LENG = {
    "5_8S_rRNA": 156, "18S_rRNA": 1869, "28S_rRNA": 2912,
    }

    # Calculate MAXLEN
    MAXLEN = int(1.2 * max(LENG.values()))

    evalue = 1E-6

    if inputfasta.endswith('.gz'):
        command = f"bgzip -@ {PBS_NCPUS} -d -k -c {inputfasta} > {PBS_JOBFS}/{sampleid}.fa"
        subprocess.run(command, check=True, shell=True)

        command = f"nhmmer --cpu {PBS_NCPUS} -E {evalue} --w_length {MAXLEN} -o /dev/null --tblout {sampleid}.primary.txt {hmmdb} {PBS_JOBFS}/{sampleid}.fa"

        print(f"Running command: {command}")
        
        # Execute nhmmer command
        subprocess.run(command, check=True, shell=True)

        os.remove(f"{PBS_JOBFS}/{sampleid}.fa")
        print(f"Removed ungzipped file for {sampleid}")

    else:
        command = f"nhmmer --cpu {PBS_NCPUS} -E {evalue} --w_length {MAXLEN} -o /dev/null --tblout {sampleid}.primary.txt {hmmdb} {inputfasta}"

        print(f"Running command: {command}")
        
        # Execute barrnap command
        subprocess.run(command, check=True, shell=True)


def read_hmm_output(hmm_output):
    #Reads hmm output into pandas df.
    # Redefine the column names for the first 12 columns
    column_names = [
        "seqid", "Target accession", "Type", "Query accession", "HMMbegin", "HMMend", "Alignstart", "Alignend", "Envstart","Envend", "Target length", "Strand", "Evalue", "Score", "Bias", "Description"
    ]

    # Reading the file with more control over parsing
    data = []
    with open(hmm_output, 'r') as file:
        for line in file:
            if not line.startswith('#'):
            # Split the line by whitespace (handles multiple spaces/tabs)
                parts = line.strip().split()
                data.append(parts)
    # Convert the list of lists to a DataFrame
    data_df = pd.DataFrame(data)
    # Assign column names to the first 12 columns
    data_df.columns = column_names + [f'Var_col_{i}' for i in range(1, len(data_df.columns) - 16 + 1)]
    # Convert relevant columns to numeric types
    numeric_columns = ["Alignstart", "Alignend", "Envstart", "Envend",
                       "HMMbegin", "HMMend", "Score", "Evalue", "Target length"]
    for col in numeric_columns:
        data_df[col] = pd.to_numeric(data_df[col])

    
    return data_df

def create_filtered_gff(sampleid, data_df):
    #Creates filtered gff file from filtered df. 
    filtered_gff = f"{sampleid}.filtered.gff"

    #Create filtered gff file
    data_gff = pd.DataFrame()

    data_gff['seqid'] = data_df["seqid"]
    data_gff["source"] = "nhmmer"
    data_gff["type"] = data_df["Type"]
    data_gff["start"] = data_df["Envstart"]
    data_gff["end"] = data_df["Envend"]
    data_gff["score"] = data_df["Evalue"]
    data_gff["strand"] = data_df["Strand"]
    data_gff["frame"] = "-"
    data_gff["attributes"] = data_df["% HMM"]

    with open(filtered_gff, 'w') as f:
        f.write("##gff-version 3\n")
        
        # Write the rows of the DataFrame to the file in GFF format
        for _, row in data_gff.iterrows():
            f.write("\t".join([str(row[col]) for col in data_gff.columns]) + "\n")

    print(f"GFF file saved to {filtered_gff}")


def process_hmm_alignments(data_df, sampleid, log_file):
    #Filters hmm alignments based on length. 
    LENG = {
    "5_8S_rRNA": 156, "18S_rRNA": 1869, "28S_rRNA": 2912,
    }

    
    data_df[['Envstart', 'Envend']] = data_df.apply(
    lambda row: (row['Envend'], row['Envstart']) if row['Strand'] == '-' else (row['Envstart'], row['Envend']),
    axis=1, result_type="expand"
    )

    data_df[['Alignstart', 'Alignend']] = data_df.apply(
    lambda row: (row['Alignend'], row['Alignstart']) if row['Strand'] == '-' else (row['Alignstart'], row['Alignend']),
    axis=1, result_type="expand"
    )

    data_df["% HMM"] = (data_df["HMMend"] - data_df["HMMbegin"] + 1) / data_df["Type"].map(LENG)
    data_df["% Target"] = (data_df["Envend"] - data_df["Envstart"] + 1) / data_df["Type"].map(LENG)


    for type in ("18S_rRNA", "5_8S_rRNA", "28S_rRNA"):
        nalignments = int(data_df[data_df["Type"] == type].shape[0])
        if nalignments == 0:
            print(f"No primary {type} for {sampleid}.Exiting")
            update_log("Errors", "Primary alignments", f"No primary {type} alignments for {sampleid}.", log_file)
            exit()
        else:
            update_log(f"rDNA_details", f"Number of primary {type} alignments", nalignments, log_file)

    #Update log file
    update_log("rDNA_details", "Number of primary alignments", int(data_df.shape[0]), log_file)
    update_log("rDNA_details", "Median Hmm % match", (data_df["% HMM"]).median(), log_file)
    update_log("rDNA_details", "Median length % match", (data_df["% Target"]).median(), log_file)

    # Filter out short alignments, retaining alignments that cover at least 95% of the target sequence
    filtered_data_df = data_df[data_df["% Target"].astype(float) >= 0.85]
    filtered_data_df = data_df[data_df["% HMM"].astype(float) >= 0.50]


    for type in ("18S_rRNA", "5_8S_rRNA", "28S_rRNA"):
        nalignments = int(filtered_data_df[filtered_data_df["Type"] == type].shape[0])
        if nalignments == 0:
            print(f"No {type} alignments pass filtering for {sampleid}.Exiting")
            update_log("Errors", "Alignment filtering", f"No {type} alignments pass filtering for {sampleid}.", log_file)
            exit()
        else:
            update_log(f"rDNA_details", f"Number of filtered {type} alignments", nalignments, log_file)


    #Update log file
    update_log("rDNA_details", "Number of filtered alignments", int(filtered_data_df.shape[0]), log_file)
    update_log("rDNA_details", "Median HMM % match for filtered alignments", (filtered_data_df["% HMM"]).median(), log_file)
    update_log("rDNA_details", "Median length % match for filtered alignments", (filtered_data_df["% Target"]).median(), log_file)

    create_filtered_gff(sampleid, filtered_data_df)

    return filtered_data_df


def morph_identification(filtered_data_df, log_file, sampleid):
    #Morph identification from filtered alignments (requires high quality 18S-28S-18S alignments)
    #Restrict to 18 and 28S for morph identification
    filtered_data_df = filtered_data_df[filtered_data_df["Type"].str.contains("18S|28S", regex=True)]

    morphs = set()
    for query_name, group in filtered_data_df.groupby("seqid"):
        group = group.sort_values(by=["Envstart"]).reset_index(drop=True)
        valid_group = set()
        i = 0
        while i < len(group):
            # start with a 18S in the plus strand, and 100bp available before this 18S on the contig
            if '18S' in group.at[i, 'Type'] and group.at[i, 'Strand'] == "+" and group.at[i, 'Envstart'] - 100 > 0:
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '28S' in group.at[i + 1, 'Type'] and group.at[i + 1, 'Strand'] == "+":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '18S' in group.at[i + 2, 'Type'] and group.at[i + 2, 'Strand'] == "+":
                        morph = group.at[i, 'seqid'], group.at[i, 'Envstart'] - 100, group.at[i + 2, 'Envstart'] - 100, '+'
                        morphs.add(morph)
            # start with a 18S in the minus strand
            elif '18S' in group.at[i, 'Type'] and group.at[i, 'Strand'] == "-":
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '28S' in group.at[i + 1, 'Type'] and group.at[i, 'Strand'] == "-":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '18S' in group.at[i + 2, 'Type'] and group.at[i + 2, 'Strand'] == "-" and group.at[i + 2, 'Envend'] + 100 <= group.at[i, 'Target length']: #Note previously we also filtered to ensure the length + 100 was < total length of query but this is not part of barrnap output
                        morph = group.at[i, 'seqid'], group.at[i, 'Envend'] + 100, group.at[i + 2, 'Envend'] + 100, '-'
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

def extract_sequences(fasta_file, morphs_df, output_dir, sampleid, log_file):
    #Extracts and saves morph sequences
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
                extracted_sequences.append((f"{seq_id}:{start}-{end}:{strand}", sequence_str))
                valid_rows.append(row)

    if len(valid_rows)==0:
        print(f"No non-ambiguous morphs for {sampleid}.Exiting")
        update_log("Errors", "Non-ambiguous morph identification", f"No morphs are non-ambiguous for {sampleid}.", log_file)
        exit()
    else:
        update_log("rDNA_details", "Number of non-ambiguous morphs", len(valid_rows), log_file)
        with open(os.path.join(output_dir, sampleid + ".rDNA.morphs.fasta"), 'w') as output_handle:
            for seq_id, sequence in extracted_sequences:
                output_handle.write(f">{seq_id}\n{sequence}\n")

        valid_morphs_df = pd.DataFrame(valid_rows)
        output_file = os.path.join(output_dir, sampleid + ".rDNA.morphs.tsv")
        valid_morphs_df.to_csv(output_file, sep='\t', index=False, header=False)

def rna_builder(sampleid, inputfasta, outputdir, log_file):
    primary_hmm = f"{sampleid}.primary.txt"

    data_df = read_hmm_output(primary_hmm)
    filtered_data_df = process_hmm_alignments(data_df, sampleid, log_file)
    morphs_df = morph_identification(filtered_data_df, log_file, sampleid)
    extract_sequences(inputfasta, morphs_df, outputdir, sampleid, log_file)

    return filtered_data_df

def get_median_morph(sampleid, PBS_NCPUS, inputfasta, log_file, filtered_data_df):
    #Reads morph sequences and finds median length unit which is designated as the representative morph for the genome. Uses this to find unit length and reference unit structure. 
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


    #Find rRNAs in refmorph
    print(f"refrdnamorph: {refrdnamorph}")

    refseqid, refstart, refend, refstrand = re.match(r"(.*):(\d+)-(\d+):([-+])", refrdnamorph).groups()
    overlapping_rRNAs = []

    overlapping_rRNAs = filtered_data_df[
        (filtered_data_df["seqid"] == refseqid) &
        (filtered_data_df["Envstart"] >= int(refstart)) &
        (filtered_data_df["Envend"] <= int(refend))
    ]

    rRNA_file = f"{sampleid}.refmorph.structure.tsv"
    with open(rRNA_file, 'w') as f:
        f.write("\t".join(filtered_data_df.columns) + "\n")
        for _, row in overlapping_rRNAs.iterrows():
            f.write("\t".join([str(row[col]) for col in overlapping_rRNAs.columns]) + "\n")

    #Use refmorph to find rDNA arrays in input genome. 
    #command=f"minimap2 -t {PBS_NCPUS} --secondary=no -o {sampleid}.asm2refmorph.paf {sampleid}.rDNA.refmorph.fasta {inputfasta}"
    print(f"Running command: {command}")

    #subprocess.run(command, shell=True, check=True)

    

def process_paf(sampleid, outputdir, log_file):
    #Process rDNA arrays to define final arrays in bed format. 
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

    result = filtered_data.groupby(["Query sequence name", "Query sequence length", "group"]).agg({"Query start": "min", "Query end": "max"}).reset_index()

    result = result[["Query sequence name", "Query sequence length", "Query start", "Query end"]]

    result = result[(result["Query end"] - result["Query start"]) > minregionlength]

    total_length = (result['Query end'] - result['Query start']).sum()

    contigs = len(pd.unique(result["Query sequence name"]))

    contigs_sorted = result.sort_values(by='Query sequence length')

    median_contig = len(contigs_sorted) // 2
    median_contig_length = contigs_sorted.iloc[median_contig]['Query sequence length']

    #Update log file
    update_log("rDNA_details", "Total length", int(total_length), log_file)
    update_log("rDNA_details", "Number of contigs", int(contigs), log_file)
    update_log("rDNA_details", "Number of rDNA arrays", int(result.shape[0]), log_file)
    update_log("rDNA_details", "Median contig length", int(median_contig_length), log_file)

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
    parser.add_argument('-l', '--hmmdb', required=True, help='Location of HMM db')
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
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fna.gz")+ 
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fa")
    )[0]
    
    run_hmm(args.ncpus, args.sampleid, args.output_dir, inputfasta, log_file, args.hmmdb)

    data_df = rna_builder(args.sampleid, inputfasta, args.output_dir, log_file)

    get_median_morph(args.sampleid, args.ncpus, inputfasta, log_file, data_df)

    process_paf(args.sampleid, args.output_dir, log_file)


 

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
