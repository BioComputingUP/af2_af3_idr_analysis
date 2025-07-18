import argparse
import os 
import sys

import pandas as pd

import utils

"""
This script parses the Disorder-PDB fasta file and saves the result in the output dir (data/parsed_fasta)
"""
def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta', help='Disorder-PDB fasta file')
    parser.add_argument('--output_dir', help='Output directory')
    parser.add_argument('--disprot_uniprot_mapping', help='Mapping between DisProt and Uniprot IDs')
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    return args



def process_protein(acc, values,disprot_uniprot_mapping):
    """
    Process the protein data and calculate disorder content, order content, and other metrics.
    """
    sequence = values["sequence"]
    labels = values["labels"]

    seq_len = len(sequence)
    disorder_count = labels.count('1')
    order_count = labels.count('0')
    disorder_content = disorder_count / seq_len
    order_content = order_count /seq_len
    na_count = labels.count('-') 
    na_content = na_count /seq_len
    annotated_content = (disorder_count + order_count) / seq_len
    category = "IDP" if disorder_content >= 0.95 else "PosReg" if disorder_content < 0.95 and disorder_content == annotated_content else "IDR"
    return {
        "DisProt_id": acc,
        "Uniprot_id": utils.extract_uniprot_id(acc,disprot_uniprot_mapping),
        "seq_len": seq_len,
        "disorder_count": disorder_count,
        "disorder_content": disorder_content,
        "order_count":order_count,
        "order_content": order_content,
        "na_count": na_count,
        "na_content": na_content,
        "annotated_content": annotated_content,
        "category": category,
        "sequence": sequence,
        "labels": labels
    }

def build_protein_dataframe(data,disprot_uniprot_mapping):
    """
    Build a dataframe with the protein data and calculate disorder content, order content, and other metrics.
    """
    protein_entries = []
    for acc, values in data.items():
        result = process_protein(acc, values, disprot_uniprot_mapping)
        protein_entries.append(result)
    protein_entries_df = pd.DataFrame(protein_entries)
    return protein_entries_df

def get_terminal(start,end,seq_len):
    if (start == 1) and (end == seq_len):
        return "IDP"
    elif start == 1 :
        return "N"
    elif end == seq_len:
        return "C"
    else:
        return "Internal"
    

def build_region_dataframe(data, disprot_uniprot_mapping):
    """
    Build a dataframe with the disorder regions for each protein.
    """
    entries = []
    for acc , values in data.items():
        sequence = values["sequence"]
        labels = values["labels"]
        segments = segment_label_regions(labels)
        for segment in segments:
            entries.append({
                "DisProt_id": acc,
                "Uniprot_id": utils.extract_uniprot_id(acc, disprot_uniprot_mapping),
                "start":segment[0],
                "end":segment[1],
                "length": segment[1] - segment[0] + 1,
                "terminal": get_terminal(segment[0], segment[1], len(sequence)),
                "sequence": sequence[segment[0]-1:segment[1]]})
    return pd.DataFrame(entries)
        # sys.exit()

def segment_label_regions(label_string, target='1'):
    """
    Segments the label string into regions of the target label (default is '1').
    """
    regions = []
    in_region = False
    start = 1

    for i, char in enumerate(label_string,start=1):
        if char == target:
            if not in_region:
                in_region = True
                start = i
        else:
            if in_region:
                regions.append((start, i - 1))
                in_region = False

    if in_region:
        regions.append((start, len(label_string)))

    return regions


    
def main(args):
    disprot_uniprot_mapping = pd.read_csv(args.disprot_uniprot_mapping,index_col='DisProt_ID')
    
    data = utils.read_fasta(args.input_fasta)     
    protein_entries_df = build_protein_dataframe(data,disprot_uniprot_mapping)
    protein_entries_df.to_csv(os.path.join(args.output_dir, 'disorder_pdb_proteins.csv'), index=False)
    print(protein_entries_df['category'].value_counts())

    disorder_regions_df = build_region_dataframe(data,disprot_uniprot_mapping)
    disorder_regions_df.to_csv(os.path.join(args.output_dir, 'disorder_pdb_regions.csv'), index=False)

if __name__ == '__main__':
    args = parser()
    main(args)