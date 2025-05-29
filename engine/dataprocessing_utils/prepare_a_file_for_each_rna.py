import sys
sys.path.extend(['../..', '.'])


from engine.training_utils.glove import get_armnet, preprocess_with_shape_armnet
from engine.R3FOLDp import ParseOneAtom

from rna_holder import RNAHolder
from pathlib import Path

import os
import pickle

import random

import torch
import torch.nn.functional as F


def prepare_rf_dataset(
        rf_pdbs_path: str | Path,
        dataset_path: str | Path) -> None:
    
    random.seed(42)
    
    pdb_filenames = sorted(os.listdir(rf_pdbs_path))
    random.shuffle(pdb_filenames)
    
    for idx, filename in enumerate(pdb_filenames):
        try:
            rnaho = RNAHolder(os.path.join(rf_pdbs_path, filename))
            ohe, pairwise_dist_matrix = rnaho.get_representation_and_eucl_matrix()
            # seq, _, _, _, _, pairwise_dist_matrix, _ = ParseOneAtom(os.path.join(rf_pdbs_path, filename))[0]
            # pairwise_dist_matrix = torch.Tensor(pairwise_dist_matrix)
            
        except Exception as err:
            print(f'Dammit! ({err})')
            continue
        
        with open(os.path.join(
                dataset_path, f'{idx + 1}.pkl'), 'wb') as file:
                pickle.dump((rnaho.sequence, pairwise_dist_matrix), file)

        if idx % 1000 == 0:
            print(f'Processedd {idx} files...')


if __name__ == '__main__':
    prepare_rf_dataset(
        rf_pdbs_path='data/rf_diffusion_data',
        dataset_path='data/rf_simple_eucl_non_baul_dataset',
    )