#!/bin/bash


# This script is used to run parse fasta file and download AF2 files

FASTA=$1
OUT_DIR=$2
DISPROT_MAPPING=$3

if [ -z "$FASTA" ] || [ -z "$OUT_DIR" ] || [ -z "$DISPROT_MAPPING" ]; then
    echo "Usage: $0 <input_fasta> <output_dir> <disprot_uniprot_mapping>"
    exit 1
fi

echo "Parsing $FASTA, saving to $OUT_DIR"

python3 parse_fasta.py --input_fasta "$FASTA" --output_dir "${OUT_DIR}/processed/parsed_fasta" --disprot_uniprot_mapping "$DISPROT_MAPPING" 
python3 download_af2_files.py --input_fasta "$FASTA" --disprot_uniprot_mapping "$DISPROT_MAPPING" --output_dir "${OUT_DIR}/AF2" --log_file ./logs/download_af2_files.log

mkdir -p "${OUT_DIR}/AF2/colab_fold"
cd "${OUT_DIR}/AF2/colab_fold"

singularity run --nv -B ~/.cache:/cache -B $(pwd):/work ../colabfold_1.5.5-cuda12.2.2.sif colabfold_batch --amber /work/afdb_failed.fasta /work/outputs/

cd outputs
find . -type f  -name "*_relaxed_rank_*.pdb" -exec cp {} "${OUT_DIR}/structures/" \; # copy the colabfold structures to directory with AF2 structures