import streamlit as st

import sys
sys.path.extend(['.', '..', '../..', 'src', '../src'])

from armnet_insides.model import ShapeArmNet, StrucArmNet
from armnet_insides.config import ArmNetConfig
from training_utils.dataflow import collate_armnetize_mask
from dataprocessing_utils.rna_holder import RNAHolder

import torch
from torch import nn

import numpy as np
import pandas as pd

from tqdm.auto import tqdm
import matplotlib.pyplot as plt
import pickle

shape_device = 'cpu'
struc_device = 'cpu'

import os
print(os.listdir('.'))

def move_dict_to_device(dct, device):
    for key in dct:
        dct[key] = dct[key].to(device)
    return dct


class CoupledModel(nn.Module):
    def __init__(self,
                 armnet_weights_path,
                 struc_model_depth: int,
                 struc_model_hidden_dim: int,
                 freeze_shape_part: bool = True,
                 shape_device: str = 'cpu',
                 struc_device: str = 'cpu') -> None:
        
        super().__init__()
        self.shape_device = shape_device
        self.struc_device = struc_device
        
        self.shape_armnet = ShapeArmNet(
            depth = ArmNetConfig.num_encoder_layers, 
             num_convs = ArmNetConfig.num_conv_layers,
             adj_ks = ArmNetConfig.conv_2d_kernel_size,
             attn_kernel_size = ArmNetConfig.conv_1d_kernel_size,
             dropout = ArmNetConfig.dropout,
             conv_use_drop1d = ArmNetConfig.conv_1d_use_dropout,
             use_bppm = ArmNetConfig.use_bppm,
        )
        self.shape_armnet.load_state_dict(torch.load(armnet_weights_path))
        self.shape_armnet.to(self.shape_device)

        if freeze_shape_part:
            for param in self.shape_armnet.parameters():
                param.requires_grad = False
        
        self.struc_armnet = StrucArmNet(depth=struc_model_depth, dim=192)
        self.struc_armnet.to(self.struc_device)

    def forward(self, tokenized_rna_batch: torch.Tensor) -> torch.Tensor:
        # print(tokenized_rna_batch)
        # tokenized_rna_batch = tokenized_rna_batch.to(self.shape_device)
        shape_representation = self.shape_armnet(tokenized_rna_batch)

        predicted_coords = self.struc_armnet(shape_representation, tokenized_rna_batch)
        # print(predicted_coords.shape)
        predicted_coords = predicted_coords[:, 1:-1, :]
        
        ## By [1:-1] we achive absence of <bos> token and the last <pad> token.
        ## <eos> token stays, but it's location is padded, so it doesn't matter
        
        return torch.cdist(predicted_coords, predicted_coords, p=2)

    def change_devices(self, 
                       shape_device: str, struc_device: str) -> 'CoupledModel':
        self.shape_device = shape_device
        self.shape_armnet.to(self.shape_device)
        
        self.struc_device = struc_device
        self.struc_armnet.to(self.struc_device)
        return self
    

class BinnedModel(nn.Module):
    def __init__(self, num_bins: int) -> None:
        super().__init__()
        self.bin_models = nn.ModuleList([nn.Module() for _ in range(num_bins)])
        self.modulo = 240 // num_bins
        self.rnaho = RNAHolder()

    def add_model(self, model_instance: nn.Module, bin_number: int) -> 'BinnedModel':
        assert bin_number < len(self.bin_models)
        self.bin_models[bin_number] = deepcopy(model_instance).cpu()
        return self

    def forward(self, x: dict, bin_number: int) -> torch.Tensor:
        assert bin_number < len(self.bin_models)
        return self.bin_models[bin_number](x)

    def change_devices(self,
                      shape_device: str | torch.device,
                      struc_device: str | torch.device,
                      bin_numbers: list[int] | None = None) -> 'BinnedModel':
        indices = range(len(self.bin_models)) if bin_numbers is None else bin_numbers
        for idx in indices:
            self.bin_models[idx].change_devices(shape_device, struc_device)

    def process_single_rna_sequence(self, seq: str, shape_device: str | torch.device = 'cuda:0') -> torch.Tensor:
        seq = seq.translate(str.maketrans('U', 'T'))
        bin_num = min(len(seq) // self.modulo, len(self.bin_models) - 1)
        
        rna_tensor_dict = collate_armnetize_mask([[seq, torch.Tensor([0]).view(1, 1, 1,)], ])[0][0]
        distance_matrix = self.forward(rna_tensor_dict, bin_num)[0].detach().cpu()

        coords = self.rnaho.calculate_coordinates_from_euclidean_distance_matrix(distance_matrix)
        fig = self.rnaho.plot_reconstr_mol_from_eucl_coords(coords)
        return fig
    

with open('finetuned_4layer_binned_model.pkl', 'rb') as file:
    binned_model = pickle.load(file)

# binned_model.process_single_rna_sequence('AAAA')

st.title('Short RNA structure prediction using light-weight transformer')

# Ввод строки пользователем
input_str = st.text_input('Enter your RNA sequence', 'AUGAGCGGCGAGCGAGCCGAGCAGCGAGC')

if input_str:
    fig = binned_model.process_single_rna_sequence(input_str)
    st.plotly_chart(fig, use_container_width=True)
