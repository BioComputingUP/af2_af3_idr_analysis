import argparse
import logging 
import os 

import pandas as pd
import requests
from tqdm import tqdm
import utils 


def download_alphafold_files(uniprot_id, out_dir=".", disprot_id=None):
    base_url = "https://alphafold.ebi.ac.uk/files"
    base_filename = f"AF-{uniprot_id}-F1"

    files = {
        "structure_cif": f"{base_filename}-model_v4.cif",
        # "pae_json": f"{base_filename}-predicted_aligned_error_v4.json"
    }
    

    for label, filename in files.items():
        url = f"{base_url}/{filename}"
        response = requests.get(url)

        if label == "structure_cif":
            output_path = f"{out_dir}/structures/{disprot_id}.cif"
        elif label == "pae_json":
            output_path = f"{out_dir}/pae/{disprot_id}.json"

        if response.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(response.content)
        else:
            logging.info(f"Failed to download {disprot_id} : {filename}: {response.status_code}")
            fails[disprot_id] = data[disprot_id]['sequence']
    

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_fasta', help='Disorder-PDB fasta file')
    parser.add_argument('--disprot_uniprot_mapping', help='Mapping between DisProt and Uniprot IDs')
    parser.add_argument('--output_dir', help='Output directory')
    parser.add_argument('--log_file', help='Log file')
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.log_file), exist_ok=True)
    os.makedirs(f"{args.output_dir}/structures", exist_ok=True)
    os.makedirs(f"{args.output_dir}/colab_fold", exist_ok=True)
    
    # os.makedirs(f"{args.output_dir}/pae", exist_ok=True)
    return args


if __name__ == "__main__":
    args = parser()
    logging.basicConfig(filename=args.log_file, level=logging.INFO)

    disprot_uniprot_mapping = pd.read_csv(args.disprot_uniprot_mapping, index_col='DisProt_ID')
    data = utils.read_fasta(args.input_fasta)

    fails = {}
    
    for disprot_id in tqdm(data.keys()):
        download_alphafold_files(utils.extract_uniprot_id(disprot_id, disprot_uniprot_mapping), args.output_dir, disprot_id)

    with open(f'{args.output_dir}/colab_fold/afdb_failed.fasta', 'w') as f:
        for fail in fails:
            f.write(f">{fail}\n{fails[fail]}\n")

    # download_alphafold_files("Q8IWJ2")