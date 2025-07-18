import argparse
import logging
import multiprocessing as mp
import os
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from vectorized_cls_metrics.vectorized_metrics import bvaluation


def parse_args():
    parser = argparse.ArgumentParser(
            prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('references_path', help='Path to the references file')

    parser.add_argument('predictions', help="directory containing prediction file(s)")

    parser.add_argument('-r', '--refList', help='reference csv file with associations between predictor and reference',
                        required=True)

    parser.add_argument('-o', '--outputDir', default='.',
                        help='directory where the output will be written (default: cwd)')

    parser.add_argument('-b', '--labels', default=None, help='filename with labels')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="DEBUG",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels >= choice will be displayed')

    args = parser.parse_args()
    return args


def set_logger(logfile, level):
    logging.basicConfig(level=level,
                        format='%(asctime)s | %(module)-13s | %(levelname)-8s | %(message)s',
                        stream=open(logfile) if logfile else sys.stderr)


def run_bvaluation(reference):
    print('**********************')
    print(reference) ## debug
    print(reference.stem) ## debug
    print('pred_references.columns: ', pred_references.columns) ## debug
    column_of_interest = [col for col in pred_references.columns if col in reference.stem][0]
    print('column_of_interest: ', column_of_interest) ## debug
    pred_for_reference = pred_references[pred_references[column_of_interest] == 1].Method.to_list()
    print('pred_for_reference: ', pred_for_reference) ## debug
    pred_paths = list(filter(lambda x: x.stem in pred_for_reference, all_pred_paths))
    print('pred_paths: ', pred_paths) ## debug
    if len(pred_paths) == 0:
        raise ValueError(f'No predictors found for reference {reference}')

    bvaluation(str(reference), pred_paths, outpath=args.outputDir, dataset=True, bootstrap=True, target=True)


excluded = ['pdb-atleast']

if __name__ == '__main__':
    args = parse_args()
    set_logger(args.log, args.logLevel)
    all_pred_paths = list(Path(args.predictions).glob('*.caid'))
    print('all_pred_paths: ', all_pred_paths) ## debug
    all_ref_paths = list(Path(args.references_path).glob('*.fasta'))
    print('all_ref_paths: ', all_ref_paths) ## debug

    # all_ref_paths = list(filter(lambda x: x.stem not in excluded, all_ref_paths))
    # print('all_ref_paths: ', all_ref_paths) ## debug

    pred_references = pd.read_csv(args.refList, sep=',')
    print(pred_references) ## debug

    os.makedirs(args.outputDir, exist_ok=True)

    with mp.Pool(processes=mp.cpu_count() - 1) as p:
        list(tqdm(p.imap(run_bvaluation, all_ref_paths), total=len(all_ref_paths), desc='Running bvaluation'))


# python3 caid.py ../data/disorder_predictions/references/disorder/ ../data/disorder_predictions/predictions/disorder/ --refList ../data/inputs/associations.csv --outputDir ../data/disorder_predictions/caid_results/

# python3 alphafold_disorder.py -i ../AF3_AF2/data/AF2/structures/ -o ../AF3_AF2/data/disorder_predictions/AF2_with_colab_fold/AlphaFold2 -f caid -dssp mkdssp
