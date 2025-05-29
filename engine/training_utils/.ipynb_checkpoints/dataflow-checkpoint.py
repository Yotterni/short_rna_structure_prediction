from engine.armnet_insides.dataset import RNA_Dataset_Pred

from pathlib import Path
from tqdm.auto import tqdm
from copy import deepcopy

import os
os.environ['ARNIEFILE'] = 'arnie_file.txt'
from arnie.bpps import bpps

import pickle
# import sys

# sys.path.append(os.path.abspath('engine/training_utils'))

import random
import numpy as np
import pandas as pd

import torch
from torch import nn
from torch.utils.data import Dataset
from torch.nn.functional import pad
from torch.utils.data import DataLoader

import typing as tp


class BinnedDataset:
    """
    Dataset of RNAs binned by their length.
    """
    def __init__(self, dataset_path: str | Path, num_bins: int, max_seq_length: int = 240) -> None:
        self.dataset_path = dataset_path
        self.num_bins = num_bins
        self.modulo = max_seq_length // self.num_bins

        self.bins = [[] for _ in range(self.num_bins)]
        
        self.translation = str.maketrans("U", "T")
        self.no_x = str.maketrans("X", "T") # try to get rid of it pls
        self.no_line = str.maketrans("-", "T")

        progress_bar = tqdm(enumerate(sorted(os.listdir(self.dataset_path))),
                            desc='Loading data into RAM...',
                            total=len(os.listdir(self.dataset_path)))

        for idx, filename in progress_bar:
            with open(os.path.join(self.dataset_path, filename), 'rb') as f:
                loaded = pickle.load(f)
                for idx, seq in enumerate(loaded[0]):
                    if loaded[1][idx].isnan().sum() > 0:
                        continue
                    
                    bin_num = min(len(seq) // self.modulo, self.num_bins - 1)
                    self.bins[bin_num].append((seq.translate(self.translation).translate(self.no_x).translate(self.no_line), loaded[1][idx]))

    def __getitem__(self, index: int,
                    ) -> tuple[torch.Tensor, tp.Any]:
        raise NotImplementedError

    def get_item(self, bin_num: int, index: int) -> tuple[str, torch.Tensor]:
        return self.bins[bin_num][index]
    
    def __len__(self) -> int:
        return sum(len(self.bins[idx]) for idx, _ in enumerate(self.bins))
    
    def __repr__(self):
        str_data = ('BinnedDataset\n' + 'Bin number  |  Amount of seqs  |  seq_len_range\n' +
            '\n'.join([f'{bin_num}  |  {len(bin)}  |  {bin_num * self.modulo, (bin_num + 1) * self.modulo}'
                       for bin_num, bin in enumerate(self.bins)]))
        return str_data
    

# def collate_and_pad(list_of_elems: list[torch.Tensor]
#             ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
#     seqs = [elem[0] for elem in list_of_elems]

#     ditstance_matrices = [elem[1] for elem in list_of_elems]
#     max_length = max(len(distance_matrix)
#                      for distance_matrix in ditstance_matrices)
    
#     armnet_attention_mask = torch.ones(len(list_of_elems), max_length)
#     loss_mask = torch.zeros(len(list_of_elems), max_length, max_length)

#     for idx, ditstance_matrix in enumerate(ditstance_matrices):
#         length = ditstance_matrix.shape[0]
#         pad_size = max_length - length

#         ditstance_matrices[idx] = pad(
#             ditstance_matrix, (0, pad_size, 0, pad_size))
        
#         armnet_attention_mask[idx, length:] = 0
#         loss_mask[idx, :length, :length] = 1
        
#     distance_batch = torch.stack(ditstance_matrices, dim=0)
#     return seqs, distance_batch, armnet_attention_mask.bool(), loss_mask

# def get_armnet(weights_path: str, device: str | torch.device) -> nn.Module:
#     config = ArmNetConfig()
#     model =  ArmNet(
#          depth = config.num_encoder_layers, 
#          num_convs = config.num_conv_layers,
#          adj_ks = config.conv_2d_kernel_size,
#          attn_kernel_size = config.conv_1d_kernel_size,
#          dropout = config.dropout,
#          conv_use_drop1d = config.conv_1d_use_dropout,
#          use_bppm = False,
#     )
#     model.load_state_dict(
#         torch.load(weights_path, map_location="cpu")
#     )
#     model = model.eval()
#     model.to(device)
#     return model


def collate_armnetize_mask(list_of_elems: list):

    ## Sequence part
    ######################################
    
    dataset = RNA_Dataset_Pred(pd.DataFrame({
        'sequence': [elem[0] for elem in list_of_elems],
        'id': list(range(len(list_of_elems))),
        # bppm_path=''
    }))

    test_dataloader = DataLoader(
        dataset = dataset, 
        batch_size=len(list_of_elems),
        num_workers=1,
        persistent_workers=True,
        shuffle=False,
        drop_last=False
    )
    
    tokenized_seqs = next(iter(test_dataloader))

    #######################################

    ## Distance matrix and loss mask part

    #######################################

    ditstance_matrices = [elem[1] for elem in list_of_elems]
    max_length = max(len(distance_matrix)
                     for distance_matrix in ditstance_matrices)
    loss_mask = torch.zeros(len(list_of_elems), max_length, max_length)

    for idx, ditstance_matrix in enumerate(ditstance_matrices):
        length = ditstance_matrix.shape[0]
        pad_size = max_length - length

        ditstance_matrices[idx] = pad(
            ditstance_matrix, (0, pad_size, 0, pad_size))
        
        loss_mask[idx, :length, :length] = 1
        
    distance_batch = torch.stack(ditstance_matrices, dim=0)

    #######################################

    ## Output

    #######################################

    return tokenized_seqs, distance_batch, loss_mask
    


class MixedDataloader:
    def __init__(self,
                 dataset: Dataset,
                 batch_size: int,
                 collate_function: 
                 tp.Callable[[list], tuple[list[str], torch.Tensor]] = 
                 collate_armnetize_mask) -> None:
        
        self.dataset = dataset
        self.batch_size = batch_size
        self.collate_function = collate_function
            
    def epoch(self, bin_number: int, first_k_elems : int | None = None) -> tp.Generator[
        tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor],
        None, None]:

        indices = list(range(len(self.dataset.bins[bin_number])))
        random.shuffle(indices)
        if first_k_elems is not None:
            indices = indices[:first_k_elems]

        for start in range(0, len(indices), self.batch_size):
            uncollated_batch = [
                self.dataset.bins[bin_number][indices[idx]] for idx in range(
                    start, min(start + self.batch_size, len(indices)))]
            yield self.collate_function(uncollated_batch)
            

    def __len__(self):
        return self.length
    
    def num_batches(self, bin_num: int):
        return len(self.dataset.bins[bin_num]) // self.batch_size


def collate_contrafoldize_mask(list_of_elems: list):

    ## Sequence part and bppms
    ######################################
    
    dataset = RNA_Dataset_Pred(pd.DataFrame({
        'sequence': [elem[0] for elem in list_of_elems],
        'id': list(range(len(list_of_elems))),
    }))

    test_dataloader = DataLoader(
        dataset = dataset, 
        batch_size=len(list_of_elems),
        num_workers=1,
        persistent_workers=True,
        shuffle=False,
        drop_last=False
    )
    
    tokenized_seqs = next(iter(test_dataloader))

    bppms = [torch.Tensor(bpps(elem[0], package='contrafold')) for elem in list_of_elems]
    bppms = torch.stack(bppms, dim=0)

    #######################################

    ## Distance matrix and loss mask part

    #######################################

    ditstance_matrices = [elem[1] for elem in list_of_elems]
    max_length = max(len(distance_matrix)
                     for distance_matrix in ditstance_matrices)
    loss_mask = torch.zeros(len(list_of_elems), max_length, max_length)

    for idx, ditstance_matrix in enumerate(ditstance_matrices):
        length = ditstance_matrix.shape[0]
        pad_size = max_length - length

        ditstance_matrices[idx] = pad(
            ditstance_matrix, (0, pad_size, 0, pad_size))
        
        loss_mask[idx, :length, :length] = 1
        
    distance_batch = torch.stack(ditstance_matrices, dim=0)

    #######################################

    ## Output

    #######################################

    return tokenized_seqs, bppms, distance_batch, loss_mask
            