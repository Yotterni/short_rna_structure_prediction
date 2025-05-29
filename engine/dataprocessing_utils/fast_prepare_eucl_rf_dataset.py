import sys
sys.path.extend(['../..', '.'])


from engine.training_utils.glove import get_armnet, preprocess_with_shape_armnet
from engine.R3FOLDp import ParseOneAtom

from rna_holder import RNAHolder
from pathlib import Path

import os
import glob
import pickle

import random

import torch
import torch.nn.functional as F

from multiprocessing import Pool

translation = str.maketrans("U", "T")
# print(txt.translate(mytable))


def prepare_rf_dataset(
        rf_pdbs_path: str | Path,
        dataset_path: str | Path,
        train_prob: float,
        val_prob: float,
        batch_size: int) -> None:
    
    random.seed(42)
    
    sequences = []
    distance_matrices = []

    fold_names = ['train', 'val', 'test']
    fold_distribution = [train_prob, val_prob, (1 - train_prob - val_prob)]
    
    # pdb_filenames = sorted(os.listdir(rf_pdbs_path))
    location = f"{rf_pdbs_path}/*.pdb"
    pdb_filenames = sorted(glob.glob(location))
    random.shuffle(pdb_filenames)
    threads = os.cpu_count()

    with Pool(threads) as pool:
        for entry in pool.imap(ParseOneAtom, pdb_filenames):
            for data in entry:
                
                seq, _, _, _, _, my, _ = data
                #RNA sequence, lowercase chars for 3'-ends (if single sequence,
                #                                           only the last char
                #                                           is lowercase)
                # for now we can put all to UPPERCASE

                seq = seq.upper().translate(translation)
                sequences.append(seq)
                # print(seq)

                #dot-bracket line, just to look at
                # print(dbn)

                #sorted list of base pairs (i,j), where i > j (0-based)
                # print(bps)

                #input matrix of pairwise distances
                # print(mx)

                #weight matrix:
                #10     - confident distances;
                # 1     - intermediate;
                #10**-6 - FW-derived shortest paths
                # print(mw)

                #output matrix of real pairwise distances
                # print(my)
                distance_matrices.append(torch.Tensor(my))

    total_num = len(sequences)
    
    train_size = int(total_num * train_prob)
    val_size = int(total_num * val_prob)
    test_size = int(total_num * (1 - train_prob - val_prob))
    
    train_slice = slice(0, train_size)
    val_slice = slice(train_size, train_size + val_size)
    test_slice = slice(train_size + val_size, train_size + val_size + test_size)
    
    with open(os.path.join(dataset_path, 'train', 'datapice.pkl', 'wb')) as file:
        pickle.dump([sequences[train_slice], distance_matrices[train_slice]], file=file)

    with open(os.path.join(dataset_path, 'val', 'datapice.pkl', 'wb')) as file:
        pickle.dump([sequences[val_slice], distance_matrices[val_slice]], file=file)

    with open(os.path.join(dataset_path, 'test', 'datapice.pkl', 'wb')) as file:
        pickle.dump([sequences[test_slice], distance_matrices[test_slice]], file=file)
                
      


if __name__ == '__main__':
    prepare_rf_dataset(
        rf_pdbs_path='data/rf_diffusion_data',
        dataset_path='data/rf_dataset_corrected_euclidean',
        train_prob=0.7,
        val_prob=0.15,
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
