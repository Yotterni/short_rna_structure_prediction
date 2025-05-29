from pathlib import Path
from tqdm.auto import tqdm
from copy import deepcopy

import os
import pickle
import sys

sys.path.append(os.path.abspath('engine/training_utils'))

import random

import torch
from torch.utils.data import Dataset
from torch.nn.functional import pad

import typing as tp

import math
import argparse
from pathlib import Path
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
import numpy as np 
import pandas as pd

from tqdm.auto import tqdm
from armnet_insides.config import ArmNetConfig
from armnet_insides.model import ArmNet
from armnet_insides.dataset import RNA_Dataset_Test
import armnet_insides

from armnet_insides.training_utils import parameter_count, seed_everything
from armnet_insides.training_utils import DeviceDataLoader


class RNADataset(Dataset):
    def __init__(self, dataset_path: str | Path) -> None:
        super().__init__()
        self.dataset_path = dataset_path
        self.files = []
        filenames = sorted(os.listdir(self.dataset_path))
        for filename in filenames:
            self.files.append(os.path.join(self.dataset_path, filename))

        self.translation = str.maketrans("U", "T")

        self.items = []
        for filename in tqdm(self.files):
            with open(filename, 'rb') as f:
                loaded = pickle.load(f)
                if len(loaded[0]) == 240:
                    seq, matrix = loaded
                    # seq = seq.translate(self.translation)
                    self.items.append([seq, matrix])
                

    def __getitem__(self, index: int) -> tuple:
        return self.items[index]

    def __len__(self) -> int:
        return len(self.items)


model =  ArmNet(
     depth = ArmNetConfig.num_encoder_layers, 
     num_convs = ArmNetConfig.num_conv_layers,
     adj_ks = ArmNetConfig.conv_2d_kernel_size,
     attn_kernel_size = ArmNetConfig.conv_1d_kernel_size,
     dropout = ArmNetConfig.dropout,
     conv_use_drop1d = ArmNetConfig.conv_1d_use_dropout,
     use_bppm = ArmNetConfig.use_bppm,
)


model.load_state_dict(
    torch.load('shape_armnet_train_test.pth',
               map_location="cpu")
)

device = 'cuda:1'
model.to(device)

def armnet_collate_fn(list_of_elems: list, device: str | torch.device = 'cuda:1'
                     ) -> tuple[torch.Tensor, torch.Tensor]:

    sequences = [elem[0] for elem in list_of_elems if len(elem[0]) == 240]
    pairwise_matrices = [elem[1] for elem in list_of_elems if elem[1].shape[0] == 240]
    pair_matrix_batch = torch.stack(pairwise_matrices, dim=0)

    df = pd.DataFrame({'sequence_id': [0] * len(sequences),
                       'sequence': sequences,
                       'id_min': [0] * len(sequences), 
                       'id_max': [len(sequences[0])] * len(sequences)})

    test_dataset = RNA_Dataset_Test(
        df,
        use_bppm = ArmNetConfig.use_bppm,
        bppm_path = None)

    test_dataloader = DeviceDataLoader(
    DataLoader(
        dataset = test_dataset, 
        batch_size = len(list_of_elems),
        num_workers = ArmNetConfig.num_workers,
        persistent_workers=True,
        shuffle=False,
        drop_last=False),
    device = device)
    
    armnet_embs = []
    with torch.no_grad():
        for x, _ in test_dataloader:
            armnet_embs.append(model(x).cpu()[:, 1:-1])
    seq_batch = torch.cat(armnet_embs, dim=0)
    return seq_batch, pair_matrix_batch

# class RAMInDataset(Dataset):
#     def __init__(self, dataset_path: str | Path,
#                  approx_num_elems_to_load: int = -1) -> None:
#         super().__init__()
#         self.dataset_path = dataset_path
#         # self.sequences = []
#         # self.objectives = []

#         self.mixed_data = {
#             '100_or_more': [],
#             'short_not_50': [],
#             '50_exactly': []
#         }

#         self.approx_num_elems_to_load = approx_num_elems_to_load

#         for idx, filename in tqdm(enumerate(sorted(os.listdir(self.dataset_path))),
#                                   desc='Loading data into RAM...',
#                                   total=len(os.listdir(self.dataset_path))):
#             with open(os.path.join(self.dataset_path, filename), 'rb') as f:
#                 loaded = pickle.load(f)

#                 for idx, seq in enumerate(loaded[0]):
#                     if len(seq) > 99:
#                         key = '100_or_more'
#                     elif len(seq) == 50:
#                         key = '50_exactly'
#                     else:
#                         key = 'short_not_50'

#                     self.mixed_data[key].append((seq, loaded[1][idx]))

#                 # self.sequences.extend(loaded[0])
#                 # self.objectives.extend(loaded[1])

#                 # if len(loaded[1]) != len(loaded[0]):
#                 #     print(len(loaded[0]), len(loaded[1]))

#             # if approx_num_elems_to_load != -1 and (idx - 1) >= approx_num_elems_to_load:
#             #     break

#     # def __getitem__(self, index: int) -> tuple[torch.Tensor, tp.Any]:
#     #     filename = self.files[index]
#     #     with open(os.path.join(self.dataset_path, filename), 'rb') as file:
#     #         loaded = pickle.load(file)
#     #     return loaded[0], loaded[1:]
        
#         # return self.preprocessed_sequences[index], self.objectives[index]

#     def __getitem__(self, index: int,
#                     ) -> tuple[torch.Tensor, tp.Any]:
#         raise NotImplementedError

#     def get_item(self, key: str, index: int) -> tuple[str, torch.Tensor]:
#         return self.mixed_data[key][index] 

#     def __len__(self) -> int:
#         return sum(len(self.mixed_data[key]) for key in self.mixed_data)


# def collate_and_pad(list_of_elems: list[torch.Tensor]
#             ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
#     seqs = [elem[0] for elem in list_of_elems]

#     ditstance_matrices = [elem[1] for elem in list_of_elems]
#     max_length = max(len(distance_matrix)
#                      for distance_matrix in ditstance_matrices)
    
#     armnet_attention_mask = torch.ones(len(list_of_elems), max_length)
#     loss_mask = torch.ones(len(list_of_elems), max_length, max_length)

#     for idx, ditstance_matrix in enumerate(ditstance_matrices):
#         length = ditstance_matrix.shape[0]
#         pad_size = max_length - length

#         ditstance_matrices[idx] = pad(
#             ditstance_matrix, (0, pad_size, 0, pad_size))
        
#         armnet_attention_mask[idx, length:] = 0
#         loss_mask[idx, length:, length:] = 0
        
#     distance_batch = torch.stack(ditstance_matrices, dim=0)
#     return seqs, distance_batch, armnet_attention_mask.bool(), loss_mask


# class MixedDataloader:
#     def __init__(self,
#                  dataset: Dataset,
#                  batch_size: int,
#                  collate_function: 
#                  tp.Callable[[list], tuple[list[str], torch.Tensor]],
#                  weights: tp.Optional[dict['str', float]] = None) -> None:
        
#         self.dataset = dataset
#         self.batch_size = batch_size
#         self.collate_function = collate_function

#         self.subdatasets = {}
#         self.length = 0

#         for subdataset_name in self.dataset.mixed_data:
#             subdataset_length = len(
#                 self.dataset.mixed_data[subdataset_name])
#             self.length += subdataset_length

#             if weights is None:
#                 subdataset_weight = subdataset_length
#             else:
#                 subdataset_weight = weights[subdataset_name]

#             self.subdatasets[subdataset_name] = [
#                 subdataset_length,
#                 subdataset_weight,
#                 list(range(subdataset_length)),
#                 0, # current index
#                 True # is not finished
#             ]
        
#         if weights is None:
#             for subdataset_name in self.subdatasets:
#                 self.subdatasets[subdataset_name][1] /= self.length
            
#     def epoch(self) -> tp.Generator[
#         tuple[torch.Tensor, torch.Tensor], None, None]:

#         epoch_mixed_indices = deepcopy(self.subdatasets)
#         for key in epoch_mixed_indices:
#             random.shuffle(epoch_mixed_indices[key][2])

#         # for key in epoch_mixed_indices:
#         #     print(epoch_mixed_indices[key][4])

#         while sum(epoch_mixed_indices[key][4] for key in epoch_mixed_indices) > 0:
#             active_subdatasets = [key for key in epoch_mixed_indices
#                                   if epoch_mixed_indices[key][4]]
            
#             current_weights = [epoch_mixed_indices[key][1] 
#                                for key in active_subdatasets]
            
#             curr_subdata = random.choices(
#                 active_subdatasets,
#                 current_weights, k=1)[0]
            
#             # print(active_subdatasets, current_weights, curr_subdata)
#             # curr_subdata = curr_subdata
            
#             curr_indices = epoch_mixed_indices[curr_subdata][2]
#             start = epoch_mixed_indices[curr_subdata][3]
#             end = start + self.batch_size

#             batch_indices = curr_indices[start : end]
#             uncollated_batch = [self.dataset.mixed_data[curr_subdata][idx]
#                                 for idx in batch_indices]
            
#             epoch_mixed_indices[curr_subdata][3] = end
#             if end >= epoch_mixed_indices[curr_subdata][0]:
#                 epoch_mixed_indices[curr_subdata][4] = False
            
#             yield self.collate_function(uncollated_batch)
            

#     def __len__(self):
#         return self.length
    
#     def num_batches(self):
#         return self.length // self.batch_size
