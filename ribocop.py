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
import shutil


def check_output(output_dir, sampleid):
    """Check if output has already been generated for a given sample ID."""
    checkpoint_file = os.path.join(output_dir, f"{sampleid}.ribocop")
    if os.path.exists(checkpoint_file + ".done"):
        print(f"Task already completed for {sampleid}. Exiting")
        sys.exit(0)
    log_file = os.path.join(output_dir, f"{sampleid}_log.json")
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            content = json.load(f)
            if 'Errors' in content:  # Check for actual key
                print(f"{sampleid} previously produced error. Exiting")
                sys.exit(0)
    open(checkpoint_file + ".running", 'w').close()



def modify_fasta(inputfasta, outputfasta, sampleid):
    #nhmmer fails when it does not detect all bases in first few hundred base pairs. can overcome this by just replacing the first 4 bp in each sequence with acgt
    with open(inputfasta, 'r') as infile, open(outputfasta, 'w') as outfile:
        for record in SeqIO.parse(infile, "fasta"):
        # Modify the sequence by adding "ACGT" to the start
            sequence = "acgt" + str(record.seq)[4:]
            # Write the modified record to the output file
            record.seq = Seq(sequence)
            SeqIO.write(record, outfile, "fasta")
        print(f"FASTA file modified for {sampleid}")

def run_hmm(PBS_NCPUS, sampleid, output_dir, inputfasta, log_file, hmmdb, tmp_dir, evalue):
    #Align the input fasta to HMM library using nhmmer.
    

    LENG = {
    "5_8S_rRNA": 154, "18S_rRNA": 1831, "28S_rRNA": 3401,
    }

    # Calculate MAXLEN
    MAXLEN = int(1.2 * max(LENG.values()))


    tmpfasta = f"{tmp_dir}/{sampleid}tmp.fa"

    if inputfasta.endswith('.gz'):
        #nhmmer requires input files to be rewindable for hmm files with multiple models (which ours has). Therefore gzipped files need to be ungzipped.
        command = f"bgzip -@ {PBS_NCPUS} -d -k -c {inputfasta} > {tmp_dir}/{sampleid}.fa"
        print(f"Running command: {command}")
        subprocess.run(command, check=True, shell=True)
        


        modify_fasta(f"{tmp_dir}/{sampleid}.fa", tmpfasta, sampleid)
        os.remove(f"{tmp_dir}/{sampleid}.fa")
        print(f"Removed ungzipped file for {sampleid}")
        
        #Run NHMMER
        command = f"nhmmer --cpu {PBS_NCPUS} -E {evalue} --w_length {MAXLEN} -o /dev/null --tblout {output_dir}/{sampleid}.primary.txt {hmmdb} {tmpfasta}"
        print(f"Running command: {command}")
        subprocess.run(command, check=True, shell=True)
        

        #remove temporary files to save space
        os.remove(f"{tmp_dir}/{sampleid}tmp.fa")
        print(f"Removed modified file for {sampleid}")

    else:
        modify_fasta(inputfasta, tmpfasta, sampleid)
        command = f"nhmmer --cpu {PBS_NCPUS} -E {evalue} --w_length {MAXLEN} -o /dev/null --tblout {output_dir}/{sampleid}.primary.txt {hmmdb} {tmpfasta}"

        print(f"Running command: {command}")
        
        # Execute barrnap command
        subprocess.run(command, check=True, shell=True)
        os.remove(f"{tmp_dir}/{sampleid}tmp.fa")
        print(f"Removed modified file for {sampleid}")


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

def create_filtered_gff(sampleid, data_df, output_dir):
    #Creates filtered gff file from filtered df. 
    filtered_gff = f"{output_dir}/{sampleid}.filtered.gff"

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


def process_hmm_alignments(data_df, sampleid, log_file, output_dir):
    #Filters hmm alignments based on length. 
    LENG = {
    "5_8S_rRNA": 154, "18S_rRNA": 1831, "28S_rRNA": 3401,
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


    for type in ("18S_rRNA", "28S_rRNA"):
        nalignments = int(data_df[data_df["Type"] == type].shape[0])
        if nalignments == 0:
            print(f"No primary {type} for {sampleid}.Exiting")
            update_log("Errors", "Primary alignments", f"No primary {type} alignments for {sampleid}.", log_file)
            exit()
        else:
            update_log(f"rDNA_details", f"Number of primary {type} alignments", nalignments, log_file)

    n_5_8S=int(data_df[data_df["Type"] == "5_8S_rRNA"].shape[0])
    update_log(f"rDNA_details", f"Number of primary 5_8S_rRNA alignments", n_5_8S, log_file)

    #Update log file
    update_log("rDNA_details", "Number of primary alignments", int(data_df.shape[0]), log_file)
    update_log("rDNA_details", "Median Hmm % match", (data_df["% HMM"]).median(), log_file)
    update_log("rDNA_details", "Median length % match", (data_df["% Target"]).median(), log_file)

    #Before filtering out short alignments, broken alignments need to be merged. 
    merged_data_df = combine_broken_alignments(data_df, sampleid, log_file)
    #recalculate %hmm percent
    merged_data_df["% HMM"] = (merged_data_df["HMMend"] - merged_data_df["HMMbegin"] + 1) / merged_data_df["Type"].map(LENG)
    merged_data_df["% Target"] = (merged_data_df["Envend"] - merged_data_df["Envstart"] + 1) / merged_data_df["Type"].map(LENG)
    #Filter out short alignments, retaining alignments that cover at least 95% of the target sequence
    filtered_data_df = merged_data_df[merged_data_df["% Target"].astype(float) >= 0.85]
    filtered_data_df = merged_data_df[merged_data_df["% HMM"].astype(float) >= 0.80]


    for type in ("18S_rRNA", "28S_rRNA"):
        nalignments = int(filtered_data_df[filtered_data_df["Type"] == type].shape[0])
        if nalignments == 0:
            print(f"No {type} alignments pass filtering for {sampleid}.Exiting")
            update_log("Errors", "Alignment filtering", f"No {type} alignments pass filtering for {sampleid}.", log_file)
            exit()
        else:
            update_log(f"rDNA_details", f"Number of filtered {type} alignments", nalignments, log_file)

    n_5_8s=int(filtered_data_df[filtered_data_df["Type"] == "5_8S_rRNA"].shape[0])
    update_log(f"rDNA_details", f"Number of filtered 5_8S_rRNA alignments", n_5_8s, log_file)


    #Update log file
    update_log("rDNA_details", "Number of filtered alignments", int(filtered_data_df.shape[0]), log_file)
    update_log("rDNA_details", "Median HMM % match for filtered alignments", (filtered_data_df["% HMM"]).median(), log_file)
    update_log("rDNA_details", "Median length % match for filtered alignments", (filtered_data_df["% Target"]).median(), log_file)

    create_filtered_gff(sampleid, filtered_data_df, output_dir)

    return filtered_data_df

def find_last_gap_envend(i, df, variable):
    if df.at[i, "type"] != "separate":
        return df.at[i, variable]

    j = i + 1
    new_value = df.at[i, variable]  # default fallback

    # Look forward for contiguous "gap" rows
    while j < len(df) and df.at[j, "type"] == "gap":
        new_value = df.at[j, variable]
        j += 1

    return new_value

def combine_broken_alignments(df, sampleid, log_file):
    #sort out -ve strand
    reset_dfs = []
    total_gaps = 0
    for combo, group in df.groupby(["seqid", "Type", "Strand"], group_keys=False):
        group = group.sort_values(by = ["Envstart"]).reset_index(drop=True)
        
        if (group["Strand"] == "+").all(): 
            group["query_gap"] = group["Envstart"] - group["Envend"].shift(fill_value=group["Envstart"].iloc[0])
            group["ref_gap"] = group["HMMbegin"] - group["HMMend"].shift(fill_value=group["HMMbegin"].iloc[0])

            conditions = [(group["ref_gap"] < -100), (group["query_gap"] > 2000), (group["ref_gap"] == 0) & (group["query_gap"] == 0)]
            choices = ["separate", "separate", "single"]

        
            group["type"] = np.select(conditions, choices, default="gap")
            total_gaps += (group["type"] == "gap").sum()
            group["merged_envend"] = [find_last_gap_envend(i, group, "Envend") for i in range(len(group))]
            group["merged_hmmend"] = [find_last_gap_envend(i, group, "HMMend") for i in range(len(group))]
            group = group[group["type"] != "gap"]
            group.loc[:, "Envend"] = group["merged_envend"]
            group.loc[:, "HMMend"] = group["merged_hmmend"]
            group = group[["seqid", "Target accession", "Type", "Query accession", "HMMbegin", "HMMend", "Alignstart", "Alignend", "Envstart","Envend", "Target length", "Strand", "Evalue", "Score", "Bias", "Description", "% HMM", "% Target"]]

            reset_dfs.append(group)
        
        else:
            group["query_gap"] = group["Envstart"] - group["Envend"].shift(fill_value=group["Envstart"].iloc[0])
            group["ref_gap"] = group["HMMbegin"].shift(fill_value=group["HMMbegin"].iloc[0]) - group["HMMend"]

            conditions = [(group["ref_gap"] < -100), (group["query_gap"] > 2000), (group["ref_gap"] == 0) & (group["query_gap"] == 0)]
            choices = ["separate", "separate", "single"]

        
            group["type"] = np.select(conditions, choices, default="gap")
            total_gaps += (group["type"] == "gap").sum()
            group["merged_envend"] = [find_last_gap_envend(i, group, "Envend") for i in range(len(group))]
            group["merged_hmmstart"] = [find_last_gap_envend(i, group, "HMMbegin") for i in range(len(group))]
            group = group[group["type"] != "gap"]
            group.loc[:, "Envend"] = group["merged_envend"]
            group.loc[:, "HMMbegin"] = group["merged_hmmstart"]
            group = group[["seqid", "Target accession", "Type", "Query accession", "HMMbegin", "HMMend", "Alignstart", "Alignend", "Envstart","Envend", "Target length", "Strand", "Evalue", "Score", "Bias", "Description", "% HMM", "% Target"]]

            reset_dfs.append(group)
    reset_df = pd.concat(reset_dfs, ignore_index=True)
    print(f"{total_gaps} merged for {sampleid}")
    update_log("rDNA_details", "Number of broken alignments", int(total_gaps), log_file)
    
    return reset_df

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
        morphs_df["morph_number"] = [f"morph{i+1}" for i in range(len(morphs_df))]
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
        morph_id = row["morph_number"]
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
                extracted_sequences.append((f"{morph_id}.{seq_id}:{start}-{end}:{strand}", sequence_str))
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

def rna_builder(sampleid, inputfasta, output_dir, log_file):
    primary_hmm = f"{output_dir}/{sampleid}.primary.txt"

    data_df = read_hmm_output(primary_hmm)
    filtered_data_df = process_hmm_alignments(data_df, sampleid, log_file, output_dir)
    morphs_df = morph_identification(filtered_data_df, log_file, sampleid)
    extract_sequences(inputfasta, morphs_df, output_dir, sampleid, log_file)

    return filtered_data_df, morphs_df

def get_median_morph(sampleid, PBS_NCPUS, inputfasta, log_file, filtered_data_df, output_dir, morphs_df):
    #Reads morph sequences and finds median length unit which is designated as the representative morph for the genome. Uses this to find unit length and reference unit structure. 
    morph_fasta = f"{output_dir}/{sampleid}.rDNA.morphs.fasta"

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

    command = f"samtools faidx {morph_fasta} {refrdnamorph} > {output_dir}/{sampleid}.rDNA.refmorph.fasta"

    print(f"Running command: {command}")

    subprocess.run(command, check=True, shell=True)


    #Find rRNAs in refmorph
    print(f"refrdnamorph: {refrdnamorph}")

    refseqid, refstart, refend, refstrand = re.match(r"(.*):(\d+)-(\d+):([-+])", refrdnamorph).groups()

    refseqid = refseqid.split("morph")[-1] 
    refseqid = refseqid.split(".", 1)[-1]

    refstart=int(refstart)
    refend=int(refend)
    overlapping_ref = []

    overlapping_ref = filtered_data_df[
        (filtered_data_df["seqid"] == refseqid) &
        (filtered_data_df["Envstart"] >= int(refstart)) &
        (filtered_data_df["Envend"] <= int(refend))
    ]

    if refstrand == "+":
        overlapping_ref["Start"] = overlapping_ref["Envstart"] - refstart
        overlapping_ref["End"] = overlapping_ref["Envend"] - refstart
    else:
        overlapping_ref["Start"] = refend - overlapping_ref["Envend"] + 1
        overlapping_ref["End"] = refend - overlapping_ref["Envstart"] + 1
    

    overlapping_ref[["Start", "End", "Type"]].to_csv(
        f"{output_dir}/{sampleid}.refmorph.structure.tsv",
        sep="\t", index=False
    )

    #Overlaps in all morphs
    all_overlaps = []

    for _, morph_row in morphs_df.iterrows():
        morph_id = morph_row["morph_number"]
        seqid = morph_row["Query sequence name"]
        start = int(morph_row["Start"])
        end = int(morph_row["End"])
        strand = morph_row["Strand"]

        overlapping = filtered_data_df[
            (filtered_data_df["seqid"] == seqid) &
            (filtered_data_df["Envstart"] >= start) &
            (filtered_data_df["Envend"] <= end)
        ].copy()

        print("Overlapping features found:", overlapping_ref.shape[0])

        if overlapping.empty:
            continue

    # Adjust rRNA coordinates relative to morph
        if strand == "+":
            overlapping["Start"] = overlapping["Envstart"] - start
            overlapping["End"] = overlapping["Envend"] - start
        else:
            overlapping["Start"] = end - overlapping["Envend"] + 1
            overlapping["End"] = end - overlapping["Envstart"] + 1

        overlapping["morph_number"] = morph_id
        
        overlapping = overlapping[[
            "morph_number",
            "Start", 
            "End", 
            "Type"
        ]]

        all_overlaps.append(overlapping)
        


    if all_overlaps:
        combined = pd.concat(all_overlaps, ignore_index=True)
        output_path = os.path.join(output_dir, f"{sampleid}.structure.gff")
        combined.to_csv(output_path, sep="\t", index=False)
        print(f"Written morph structure file with {len(combined)} rRNAs from {len(morphs_df)} morphs.")
    else:
        print(f"No overlapping rRNAs found for any morph in sample {sampleid}")

    #Use refmorph to find rDNA arrays in input genome. 
    command=f"minimap2 -t {PBS_NCPUS} --secondary=no -o {output_dir}/{sampleid}.asm2refmorph.paf {output_dir}/{sampleid}.rDNA.refmorph.fasta {inputfasta}"
    print(f"Running command: {command}")

    subprocess.run(command, shell=True, check=True)

    

def process_paf(sampleid, output_dir, log_file):
    #Process rDNA arrays to define final arrays in bed format. 
    column_names = [
        "Query sequence name", "Query sequence length", "Query start", "Query end", "Relative strand",
        "Target sequence name", "Target sequence length", "Target start", "Target end", "Number of residue matches",
        "Alignment block length", "Mapping quality"
    ]

    #Read in paf file as pandas df and set output file
    data = pd.read_csv(f"{output_dir}/{sampleid}.asm2refmorph.paf", sep='\t', usecols=list(range(0,12)), names=column_names)
    outputfile = f"{output_dir}/{sampleid}.rdnaregions.bed"

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
    parser.add_argument('-p', '--species', required=False, help='Species name')
    parser.add_argument('-i', '--input_dir', required=False, default=".", help="Input directory (default: current directory)")
    parser.add_argument('-if', '--input_fasta', required=False, help="Full path to input FASTA file")
    parser.add_argument('-d', '--tmp_dir', required=False, help="Temporary directory for storing modified FASTA files (default:output directory)")
    parser.add_argument('-e', '--e_value', required=False, default=1E-6, help="E value threshold for nhmmer")

    args = parser.parse_args()



    tmp_dir = args.tmp_dir if args.tmp_dir else args.output_dir

#Create output directory if it does not exist.
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"Directory '{args.output_dir}' created successfully.")
    except OSError as e:
        print(f"Error creating directory '{args.output_dir}': {e}")
        sys.exit(1)

#Check if output already exists
    check_output(args.output_dir, args.sampleid)


    log_file = os.path.join(args.output_dir, f"{args.sampleid}_log.json")

    sampleinfo = {
        "Sample_details":{
        "Sample ID": args.sampleid,
        **({"Species": args.species} if args.species else {})
        }
    }

    with open(log_file, "w") as f:
        json.dump(sampleinfo, f, indent=2)

    if args.input_fasta:
        inputfasta = args.input_fasta
        if not os.path.exists(inputfasta):
            sys.exit(f"Error: Specified input FASTA file does not exist: {inputfasta}")
    else:
        fasta_candidates = (
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fna") + 
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fna.gz")+ 
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fa") + 
        glob.glob(f"{args.input_dir}/*/{args.sampleid}*.fasta")
        )
        if not fasta_candidates:
            sys.exit(f"Error: No FASTA file found matching sample ID {args.sampleid} in {args.input_dir}")

        inputfasta = fasta_candidates[0]
    
    run_hmm(args.ncpus, args.sampleid, args.output_dir, inputfasta, log_file, args.hmmdb, tmp_dir, args.e_value)

    filtered_data_df, morphs_df = rna_builder(args.sampleid, inputfasta, args.output_dir, log_file)

    get_median_morph(args.sampleid, args.ncpus, inputfasta, log_file, filtered_data_df, args.output_dir, morphs_df)

    process_paf(args.sampleid, args.output_dir, log_file)

    checkpoint_file = os.path.join(args.output_dir, f"{args.sampleid}.ribocop")
    
    open(checkpoint_file + ".done", 'w').close()

    os.remove(f"{checkpoint_file}.running")


    print("Task completed. Exiting")
    
    

if __name__ == "__main__":
    main()     
