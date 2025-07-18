# Modelling Intrinsically Disordered Regions from AlphaFold2 to AlphaFold3

This repository contains code for comparing the performance of AlphaFold3 and AlphaFold2 on intrinsic disorder prediction using the CAID3's Disorder-PDB dataset. 

## Data preparation

All the datapreparation including reading the fasta file, parsing the data, downloading AF2 stuctures from AFDB or generating colabfold structures could be done by running: 

    cd src
    ./prepare_data_pipeline.sh \
         <input_fasta> # Input fasta file with [disorder] labels
         <output_dir> # The directory for writing the results from fasta file, downloading/generating AF2/AF3, CAID IDR evaluations, saving figures from analysis. we used path/to/project/data as output_dir
         <disprot_uniprot_mapping>" # a csv file containing DisProt ids and their Corresponding UniProt ids. This file was obtained from DisProt Curators, or it could be extracted from DisProt directly for proteins that are released.

Individual steps could be performed as below:


### 1. Parse the fasta file
The fasta file contains protein ids, protein sequence, and the disorder label. The dataset used in this paper is the Disorder-PDB dataset from CAID3 challenge (https://caid.idpcentral.org/challenge/results). 

Format of the file should be: 


    >DP0123
    MKTFFVLLLCTFTVLSSGLTQGAE
    111111000000------111100

where 1 indicates disorder label, 0 indicates order label, and - indicates residues that are not annotated and are excluded from assessment. For CAID3 input and output formats, see https://caid.idpcentral.org/challenge.  

    python3 parse_fasta.py --input_fasta <input.fasta> --output_dir <output_dir> --disprot_uniprot_mapping <disprot_uniprot_mapping.csv>

### 2. Obtain AlphaFold2 structures

The AlphaFold2 structures used for this paper is avaialable at `protein.bio.unipd.it/shared/` and could be obtained by: 

    mkdir -p <output_dir>/AF2/
    cd <output_dir>/AF2/
    wget https://protein.bio.unipd.it/shared/caid3_disorder_pdb_af2.zip -O temp.zip && unzip temp.zip && rm temp.zip


If you want to obtain them for your own fasta file, you could run these codes:  

    # downloads AF2 stuctures from AlphaFoldDB and writes the fails in a fasta file

    python3 download_af2_files.py --input_fasta <input.fasta> --disprot_uniprot_mapping <disprot_uniprot_mapping.csv> --output_dir <output_dir> --log_file <log_dir/download_af2_files.log>


To generate ColabFold structures for structures that were not found in AlphaFoldDB, we ran colabfold in singularity container. The installation guide is provided [here](https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker). 

    cd <output_dir>
    mkdir -p AF2/colabfold
    cd AF2/colabfold

    singularity run --nv -B ~/.cache:/cache -B $(pwd):/work ../colabfold_1.5.5-cuda12.2.2.sif colabfold_batch --amber /work/afdb_failed.fasta /work/outputs/

    cd outputs
    find . -type f  -name "*_relaxed_rank_*.pdb" -exec cp {} <output_dir>/AF2/structures/ \; # copy the colabfold structures to directory with AF2 structures
    


## Obtain AlphaFold3 structures
AlphaFold3 stuctures were manually downloaded from https://alphafoldserver.com for sequences in Disrder-PDB dataset. The mmCIF files were saved in `<output_dir>/AF3/structures`. We are sharing the AlphaFold3 structures that we downloaded for Disorder-PDB dataset available at `protein.bio.unipd.it/shared/` in accordance with [AlphaFold Server Terms of Service](https://alphafoldserver.com/terms), subjected to to [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms). 

    mkdir -p <output_dir>/AF3/
    cd <output_dir>/AF3/

    wget https://protein.bio.unipd.it/shared/caid3_disorder_pdb_af3.zip -O temp.zip && unzip temp.zip && rm temp.zip
    
## Running AlphaFold-disorder package

AlphaFold disorder package could be cloned separately from [AlphaFold-disorder](https://github.com/BioComputingUP/AlphaFold-disorder/tree/main), The package needs dssp installed: https://github.com/PDB-REDO/dssp. If you have problems with mmcif_ma, try these where you clone AlphaFold-disorder:

    wget https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_ma.dic -O mmcif_ma.dic
    
    os.environ["MMCIF_MA_DIC"] = "/path/to/AlphaFold-disorder/mmcif_ma.dic"

    sudo curl -o /var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
    sudo curl -o /var/cache/libcifpp/mmcif_ma.dic https://github.com/ihmwg/ModelCIF/raw/master/dist/mmcif_ma.dic
    sudo curl -o /var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic

if you have a problem with cmake: 

    pip install --upgrade cmake
    /path/to/cmake -S . -B build
    sudo /path/to/cmake --build build
    sudo /path/to/cmake --install build

    sudo curl -o /var/cache/libcifpp/components.cif https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
    sudo curl -o /var/cache/libcifpp/mmcif_ma.dic https://github.com/ihmwg/ModelCIF/raw/master/dist/mmcif_ma.dic
    sudo curl -o /var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic


AlphaFold-disorder package generates 4 outputs: 
1. Disorder prediction for given structures, based on 1 - pLDDT
2. Disorder prediction for given structures, based on window averaged RSA
3. Binding prediction for given structures, based on combining pLDDT and RSA
4. A tsv (intermediate) file, containing `aa, plddt, 1-plddt, rsa, ss`, where ss is the secondary structure generated from dssp

To run the AlphaFold-disroder package on AF2 and AF3 structures: 

    python3 alphafold_disorder.py -i <output_dir>/AF2/structures -o <output_dir>/disorder_prediction/AF2/AlphaFold2.tsv -dssp mkdssp -f caid
    python3 alphafold_disorder.py -i <output_dir>/AF3/structures -o <output_dir>/disorder_prediction/AF3/AlphaFold3.tsv -dssp mkdssp -f caid

    # Copy the outputs of AlphaFold-disorder package to another folder that is ready for assessment

    cp <output_dir>/disorder_prediction/AF2/AlphaFold2_disorder.dat <output_dir>/disorder_prediction/predictions/disorder/AlphaFold2-plddt.caid
    cp <output_dir>/disorder_prediction/AF2/AlphaFold2_disorder-25.dat <output_dir>/disorder_prediction/predictions/disorder/AlphaFold2-rsa.caid
 
    cp <output_dir>/disorder_prediction/AF3/AlphaFold3_disorder.dat <output_dir>/disorder_prediction/predictions/disorder/AlphaFold3-plddt.caid
    cp <output_dir>/disorder_prediction/AF3/AlphaFold3_disorder-25.dat <output_dir>/disorder_prediction/predictions/disorder/AlphaFold3-rsa.caid

## Evaluation of Predictions as in CAID
caid.py script runs the evaluation of prediction files as in CAID challenge, and saves the results in the folder given in --outputDir. The evaluation code could be cloned from https://github.com/marnec/vectorized_cls_metrics. We included the code in this repository for easier access too. 

for caid.py, the disorder predictions (.caid files) must be in `<output_dir>/disorder_prediction/predictions/disorder` , and the references must be in `<output_dir>/disorder_prediction/references/disorder`. references, are the same fasta files with labels that were used in data preparation step. 

caid.py also needs a --refList, that declares which methods should be assessed on what references. For example: 

    Method,disorder,binding,linker
    AlphaFold2-rsa,1,0,0
    AlphaFold2-plddt,1,0,0
    AlphaFold3-rsa,1,0,0
    AlphaFold3-plddt,1,0,0
    AlphaFold2-binding,0,1,0
    AlphaFold3-binding,0,1,0

To execute the script: 

    cd src
    python3 caid.py \
        <output_dir>/disorder_prediction/references/disorder/ \
        <output_dir>/disorder_prediction/predictions/disorder/ \
        --refList <output_dir>/inputs/associations.csv \
        --outputDir <output_dir>/disorder_prediction/caid_results/

## Analysis
To generate paper's figures and analysis, you can run the src/figures.ipynb, src/caid_evaluations.ipynb
