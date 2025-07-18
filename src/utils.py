import argparse
import os
import sys
import subprocess
import json 
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB import MMCIFParser, PDBParser
import numpy as np
import pandas as pd
import utils
import re


def read_fasta(fasta_file):
    """
    Reads a fasta file (Disorder-PDB) and returns a dictionary with the sequence and labels for each entry.
    """
    with open(fasta_file, 'r') as f:
        lines = f.read().splitlines()

        data = {}
        current_id = ""
        sequence_lines = []
        label_line = ""

        for line in lines:
            if line.startswith(">"):
                if current_id:
                    # Save the previous entry
                    data[current_id] = {
                        "sequence": ''.join(sequence_lines),
                        "labels": label_line
                    }
                # Start a new entry
                current_id = line[1:].strip()
                sequence_lines = []
                label_line = ""
            elif set(line.strip()) <= {"0", "1", "-"}:
                label_line = line.strip()
            else:
                sequence_lines.append(line.strip())

        if current_id:
            data[current_id] = {
                "sequence": ''.join(sequence_lines),
                "labels": label_line
            }

    return data
    
def extract_uniprot_id(acc, disprot_uniprot_mapping):
    """
    Extracts the Uniprot ID from the DisProt ID.
    """
    # Check if the DisProt ID exists in the mapping
    if acc in disprot_uniprot_mapping.index:
        return disprot_uniprot_mapping.loc[acc]['UniProt_ID']
    else:
        print(f"Warning: {acc} not found in mapping.")
        return None
    
    

def get_ca_coords(structure):
    model = structure[0]
    coords = [atom.get_coord() for atom in model['A'].get_atoms() if atom.get_name() == 'CA']
    return np.array(coords)
    
def superimpose(segment_1, segment_2):
    """
    return the new coordinates of segment_2
    upon superposition with segment_1
    """
    sup = SVDSuperimposer()
    sup.set(segment_1, segment_2)
    sup.run()
    rot, tran = sup.get_rotran()
    return np.dot(segment_2, rot) + tran 


def get_parser(file_name):
    if '.pdb' in file_name:
        parser = PDBParser(QUIET=True)
    else:
        parser = MMCIFParser(QUIET=True)
    return parser

def get_distance_per_residue(acc, file1, file2 ):
    structure1 = get_parser(file1).get_structure(acc, file1)
    structure2 = get_parser(file2).get_structure(acc, file2) 
    ca_atoms_1 = get_ca_coords(structure1)  
    ca_atoms_2 = get_ca_coords(structure2)
    new_ca_atoms_2 = superimpose(ca_atoms_1,ca_atoms_2)
    rmsd = np.sqrt(np.sum((ca_atoms_1 - new_ca_atoms_2)**2)/ca_atoms_1.shape[0])
    residue_distances = np.sqrt(np.sum((ca_atoms_1 - new_ca_atoms_2)**2 , axis = 1))
    return round(rmsd,3) , residue_distances

def validate_result(dictionary):
    for key ,value in dictionary.items():
        if value == '':
            return False
    return True

def run_tmscore(file1, file2):
    acc = os.path.basename(file1).split('.')[0]
    result = subprocess.run(["TMscore", file1, file2], capture_output=True, text=True)
    output = result.stdout.strip()

    results = {
        "acc": acc,
        "TM-score": "",
        "RMSD": "",
        "RMSD_folded": "",
        "GDT_TS": "",
        "alignment_length": "",
        "common_res_length": "",
        "alignment_coverage": ""
    }

    results["alignment_length"] = output.splitlines()[-3].count(':')

    for line in output.splitlines():
        if "Number of residues in common=" in line:
            results["common_res_length"] = int(line.strip().split()[-1])
        elif "TM-score    =" in line:
            results["TM-score"] = float(line.split()[2])
        elif "RMSD of  the common residues="  in line:
            results["RMSD"] = float(line.strip().split()[-1])
        elif "GDT-TS-score=" in line:
            results["GDT_TS"] = float(line.strip().split()[1])
        elif "Superposition in the TM-score:" in line:
            match = re.search(r'Length\(d<5\.0\)=\s*([*\d]+)\s*RMSD=\s*([\d.]+)', line.strip())
            if match:
            
                results["RMSD_folded"] = float(match.group(2))
        
                results["alignment_coverage"] = round(results["alignment_length"] / results["common_res_length"],3)
      
    if validate_result(results):
        return results
    else: 
        print("TMscore output couldn't be parsed correctly, the results are:")
        print(results)
    


def calculate_rg(acc,file, start = None, end = None):
    parser = get_parser(file)
    structure = parser.get_structure(acc, file)
    model = structure[0]

    coords = np.array([atom.get_coord() for atom in model['A'].get_atoms() if atom.get_name() == 'CA'])
    coords = coords[start-1:end]
    center = np.mean(coords, axis=0)
    # print(acc, file, start, end, coords.shape )

    # Squared distances from center
    sq_dists = np.sum((coords - center) ** 2, axis=1)
    rg = np.sqrt(np.mean(sq_dists))
    return rg


def load_pae(file, mode='af2'):
    """
    Load the predicted aligned error (PAE) from a JSON file.
    mode: 'af2' or 'af3'
    """
    if mode == 'af2':
        with open(file, "r") as f:
            data = json.load(f)
        data = np.array(data[0]['predicted_aligned_error'])
  
    elif mode == 'af3':
        with open(file, "r") as f:
            data = json.load(f)
        data = np.array(data['pae'])
        
    return data
    
def get_interaction_score(data,start_1,end_1,start_2,end_2):
    """
    Get the median PAE value between 2 regions.
    The interaction score is defined as the median PAE value between the two regions, used in TED paper.
    """
    region_pae_1 = data[start_1 - 1:end_1, start_2 - 1:end_2]
    region_pae_2 = data[start_2 - 1:end_2, start_1 - 1:end_1]
    region_pae = np.concatenate((region_pae_1.flatten(), region_pae_2.flatten()))
    median_pae = np.median(region_pae)
    return median_pae