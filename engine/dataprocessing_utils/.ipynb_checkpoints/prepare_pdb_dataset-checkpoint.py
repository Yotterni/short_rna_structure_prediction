from pathlib import Path

import pickle
import os

import pandas as pd
import matplotlib.pyplot as plt
import torch


def create_kaggle_dataset(path_to_kaggle_folder: str | Path,
                          path_to_kaggle_dataset: str | Path) -> None:
    train_labels = pd.read_csv(os.path.join(path_to_kaggle_folder, 'train_labels.v2.csv'))
    val_labels = pd.read_csv(os.path.join(path_to_kaggle_folder, 'validation_labels.csv'))
    
    train_tuple = extract_seqs_and_distance_matrices(train_labels)
    val_tuple = extract_seqs_and_distance_matrices(val_labels)
    
    with open(os.path.join(path_to_kaggle_dataset,
                           'train', 'kaggle_train.pkl'), 'wb') as file:
        pickle.dump(train_tuple, file)
        
    with open(os.path.join(path_to_kaggle_dataset,
                           'val', 'kaggle_val.pkl'), 'wb') as file:
        pickle.dump(val_tuple, file)

def extract_seqs_and_distance_matrices(coord_df) -> tuple[list[str], list[torch.Tensor]]:
    list_of_seqs = []
    list_of_distance_matrices = []
    
    for name, group in coord_df.groupby(coord_df['ID'].apply(lambda x: ''.join(x.split('_')[:-1]))):
        seq = ''.join(group['resname'].tolist())
        if len(seq) < 1000:
            list_of_seqs.append(seq[:-1])
            coordinates = torch.Tensor(group[['x_1', 'y_1', 'z_1']].values)[:-1]
            list_of_distance_matrices.append(torch.cdist(coordinates, coordinates))
    return list_of_seqs, list_of_distance_matrices


if __name__ == '__main__':
    create_kaggle_dataset('data/kaggle_data/', 'data/kaggle_dataset_v2/')