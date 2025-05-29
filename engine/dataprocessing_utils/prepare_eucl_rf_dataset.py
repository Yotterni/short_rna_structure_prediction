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
        dataset_path: str | Path,
        train_prob: float,
        val_prob: float,
        batch_size: int) -> None:
    
    random.seed(42)
    
    sequences = []
    # ohe_representations = []
    distance_matrices = []

    fold_names = ['train', 'val', 'test']
    fold_distribution = [train_prob, val_prob, (1 - train_prob - val_prob)]
    
    pdb_filenames = sorted(os.listdir(rf_pdbs_path))
    random.shuffle(pdb_filenames)
    
    for idx, filename in enumerate(pdb_filenames):
        try:
            # rnaho = RNAHolder(os.path.join(rf_pdbs_path, filename))
            # ohe, pairwise_dist_matrix = rnaho.get_representation_and_eucl_matrix()
            seq, _, _, _, _, pairwise_dist_matrix, _ = ParseOneAtom(os.path.join(rf_pdbs_path, filename))[0]
        except Exception as err:
            print(f'Dammit! ({err})')
            continue
            
        sequences.append(seq.upper())
        # ohe_representations.append(ohe)
        distance_matrices.append(torch.Tensor(pairwise_dist_matrix))
                
        if len(sequences) % batch_size == 0:

            fold = random.choices(fold_names, fold_distribution, k=1)[0]

            with open(os.path.join(
                dataset_path, fold, f'{idx + 1}.pkl'), 'wb') as file:
                pickle.dump((sequences,
                             # ohe_representations,
                             distance_matrices), file)

            sequences = []
            # ohe_representations = []
            distance_matrices = []

        if idx % (batch_size) == 0:
            print(f'Processedd {idx} files...')


if __name__ == '__main__':
    prepare_rf_dataset(
        rf_pdbs_path='data/rf_diffusion_data',
        dataset_path='data/rf_dataset_corrected_euclidean',
        train_prob=0.9,
        val_prob=0.05,
        batch_size=4096
    )

# model = get_armnet('shape_armnet_train_test.pth', 'cuda:3')
# model.eval()

# def run_armnet_preprocessing(sequences: list[str],
#                              model,
#                              batch_size: int = 512,
#                              armnet_device: str | torch.device) -> None:
#     preprocessed_sequences = []
#     # for start in tqdm(range(0, len(sequences), batch_size),
#     #                   desc='Running shape-armnet preprocessing'):

#     for start in range(0, len(sequences), batch_size):
#         end = start + batch_size
#         batch_sequences = sequences[start:end]
#         processed_batch_sequences = preprocess_with_shape_armnet(
#             model,
#             batch_sequences,
#             device = armnet_device,
#         )
#         preprocessed_sequences.extend(processed_batch_sequences)
#     return preprocessed_sequences


# def prepare_rf_dataset(rf_pdbs_path: str | Path,
#                        dataset_path: str | Path,
#                        armnet_weights_path: str | Path,
#                        net_batch_size: int,
#                        file_batch_size: int,
#                        armnet_device: str | torch.device) -> None:
    
#     assert file_batch_size > net_batch_size and file_batch_size % net_batch_size == 0

#     counter = 0
#     sequences = []
#     current_sequences = []
#     preprocessed_sequences = []
#     distance_matrices = []
    
#     pdb_filenames = sorted(os.listdir(rf_pdbs_path))
#     for idx, filename in enumerate(pdb_filenames):
#         rnaho = RNAHolder(os.path.join(rf_pdbs_path, filename))

#         try:
#             _, pairwise_dist_matrix = rnaho.get_representation_and_matrix()
#         except Exception as err:
#             print(f'Dammit! ({err})')
#             continue

#         if len(rnaho.sequence) != 240:
#             delta = 240 - len(rnaho.sequence)
#             rnaho.sequence += 'A' * delta
#             pairwise_dist_matrix = (
#                 F.pad(F.pad(pairwise_dist_matrix, (0, delta)).T, (0, delta)).T
#             )
            

#         current_sequences.append(rnaho.sequence)
#         sequences.append(rnaho.sequence)
#         distance_matrices.append(pairwise_dist_matrix)

#         if len(current_sequences) % net_batch_size == 0:
#             preprocessed_sequences.extend(run_armnet_preprocessing(
#                 current_sequences,
#                 armnet_weights_path,
#                 net_batch_size,
#                 armnet_device
#             ))

#             current_sequences = []
                
#         if len(sequences) % file_batch_size == 0:
#             with open(os.path.join(dataset_path, f'{counter}.pkl'), 'wb') as file:
#                 pickle.dump((sequences,
#                              preprocessed_sequences,
#                              distance_matrices), file)
#             counter += 1

#             sequences = []
#             preprocessed_sequences = []
#             distance_matrices = []


#         if idx % (net_batch_size * 1) == 0:
#             print(f'Processedd {idx} files...')

# if __name__ == '__main__':
#     prepare_rf_dataset(
#         '../../data/rf_diffusion_data',
#         '../../data/rf_dataset',
#         '../../../V2/ArmNet_part/train_test_pseudolabel_weights/model.pth',
#         256,
#         4096,
#         'cuda:1')
