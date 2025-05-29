import sys
sys.path.append('.')

from engine.dataprocessing_utils.rna_holder import RNAHolder
from engine.training_utils.dataflow import BinnedDataset, MixedDataloader, collate_and_pad
from engine.training_utils.glove import get_armnet, preprocess_with_shape_armnet

from torch.utils.data import DataLoader

from tqdm import tqdm

# import engine
import armnet_insides

import torch
import torch.nn as nn
import torch.nn.functional as F

import matplotlib.pyplot as plt
import datetime
import pickle

import warnings
warnings.filterwarnings('ignore')

struc_device = 'cuda:0'
shape_device = 'cuda:1'


rf_train = BinnedDataset('data/rf_dataset/train', 5)
rf_val = BinnedDataset('data/rf_dataset/val', 5)
rf_test = BinnedDataset('data/rf_dataset/test', 5)


train_loader = MixedDataloader(rf_train, 256, collate_function=collate_and_pad)
val_loader = MixedDataloader(rf_val, 256, collate_function=collate_and_pad)
test_loader = MixedDataloader(rf_test, 256, collate_function=collate_and_pad)

def get_loaders(batch_size: int) -> tuple:
    train_loader = MixedDataloader(rf_train, batch_size, collate_function=collate_and_pad)
    val_loader = MixedDataloader(rf_val, batch_size, collate_function=collate_and_pad)
    test_loader = MixedDataloader(rf_test, batch_size, collate_function=collate_and_pad)
    return train_loader, val_loader, test_loader

shape_armnet = get_armnet(
    'shape_armnet_train_test.pth',
    device=shape_device
)


def process_batch(
        batch: tuple[list[str], torch.Tensor, torch.Tensor]
        ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    sequences, distance_matrix_batch, attention_mask, loss_mask = batch

    shape_emb_batch = preprocess_with_shape_armnet(
        shape_armnet,
        sequences=sequences,
        device=shape_device,
    )

    shape_emb_batch = shape_emb_batch.to(struc_device)
    distance_matrix_batch = distance_matrix_batch.to(struc_device)

    distance_matrix_batch /= 100 # TODO: убрать!!!!!!!!
    attention_mask = attention_mask.to(struc_device)
    # attention_mask = None
    loss_mask = loss_mask.to(struc_device)

    # print('shape emb batch shape:', shape_emb_batch.shape)
    # print('distance_matrix_batch shape:', distance_matrix_batch.shape)
    # print('attention_mask shape:', attention_mask.shape)

    return shape_emb_batch, distance_matrix_batch, attention_mask, loss_mask 


def train_step(model: nn.Module,
          loader: DataLoader,
          bin_num: int,
          criterion: nn.Module,
          optimizer: torch.optim.Optimizer,
          scheduler: torch.optim.lr_scheduler,
          history_loss: list[float]) -> None:
    pbar = tqdm(loader.epoch(bin_num), total=loader.num_batches(bin_num), leave=True)

    for batch in pbar:
        shape_emb_batch, distance_matrix_batch, attention_mask, loss_mask = process_batch(batch)
        predicted_matrix = model(shape_emb_batch, attention_mask)
        # print(pair_matrix_batch.shape)
        loss = criterion(predicted_matrix, distance_matrix_batch, loss_mask)
        history_loss.append(loss.item())
        loss.backward()
        optimizer.step()
        optimizer.zero_grad()
        scheduler.step()
        pbar.set_postfix_str(f'loss = {loss.item():.4f}')


@torch.no_grad()
def val_step(model: nn.Module,
         loader: DataLoader, 
         bin_num: int,
         criterion: nn.Module,
         history_loss: list[float]) -> None:
    pbar = tqdm(loader.epoch(bin_num), total=loader.num_batches(bin_num), leave=True)
    for batch in pbar:

        shape_emb_batch, distance_matrix_batch, attention_mask, loss_mask = process_batch(batch)
        predicted_matrix = model(shape_emb_batch, attention_mask)
        loss = criterion(predicted_matrix, distance_matrix_batch, loss_mask)

        history_loss.append(loss.item())
        pbar.set_postfix_str(f'loss = {loss.item():.4f}')


class MaskedLoss(nn.Module):
    def __int__(self, mode: str):
        super().__init__()

        self.mode = mode

    def forward(self, x_predicted: torch.Tensor,
                x_true: torch.Tensor,
                loss_mask: torch.Tensor) -> torch.Tensor:
        """
        Expects loss mask to be 1 when this position is padded and 0 otherwise.
        """
        return ((x_predicted - x_true) ** 2 * loss_mask).sum() / loss_mask.sum()


def run_experiment(model: nn.Module,
                   train_loader: DataLoader,
                   val_loader: DataLoader,
                   num_epochs: int,
                   criterion: nn.Module,
                   optimizer: torch.optim.Optimizer,
                   scheduler: torch.optim.lr_scheduler,
                   history_train: list[float],
                   history_val: list[float]) -> tuple[nn.Module, float, float]:

    batch_sizes = [1024, 512, 256, 128, 96, 96]
    for bin_num in range(2, rf_train.num_bins):
        print(f'Bin num: {bin_num}\n')

        train_loader, val_loader, _ = get_loaders(batch_sizes[bin_num])

        for epoch in range(num_epochs):
            print(f'Epoch: {epoch}')
            print(datetime.datetime.now())
            train_step(model, train_loader, bin_num, criterion, optimizer, scheduler, history_train)
            val_step(model, val_loader, bin_num, criterion, history_val)
            print(f'Train loss: {history_train[-1]:.04f}')
            print(f'Val loss: {history_val[-1]:.04f}')
            with open(f'results_for_bin_{bin_num}.pkl', 'wb') as file:
                    pickle.dump([model, history_train_loss, history_val_loss], file=file)

    print(f'Train loss: {history_train[-1]:.04f}')
    print(f'Val loss: {history_val[-1]:.04f}')

    # nucleotides, matrix, attn_mask = next(iter(train_loader))
    # model.plot_predicted_structure(model(nucleotides.to(device), attn_mask)[0].detach().cpu())
    # model.plot_predicted_structure(matrix[0])
    return model, history_train, history_val


class ArmnetLikeModel(armnet_insides.model.ArmNet):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.rnaho = RNAHolder()
        self.initial_projection = nn.Sequential(
            nn.Linear(192, kwargs['dim']),
            nn.ReLU()) if 'dim' in kwargs else nn.Identity()
            
        # self.final_projection_list = nn.ModuleList([
        #     nn.Linear(kwargs['dim'], kwargs['dim'] * 4),
        #     nn.ReLU(),
        #     nn.LayerNorm(kwargs['dim'] * 4),
        #     nn.Linear(kwargs['dim'] * 4, kwargs['dim'] * 4),
        #     nn.ReLU(),
        #     nn.LayerNorm(kwargs['dim'] * 4),
        #     nn.Linear(kwargs['dim'] * 4, 5)]) # TODO: тюнить эту чиселку

        self.final_projection = nn.Linear(kwargs['dim'], 4)

    def forward(self, x: torch.Tensor, attn_mask: torch.Tensor) -> torch.Tensor:
        x = self.initial_projection(x)
        # mask = x['forward_mask']
        bpp_mask = None
        adj = self.outer_product_mean(x)
        # print(attn_mask.shape)
        x = self.transformer(x, adj, mask=attn_mask, bpp_mask=bpp_mask) # 

        x = self.final_projection(x)
        # return torch.cdist(x, x, p=2).view(x.shape[0],
        #                                    x.shape[1],
        #                                    x.shape[1])
        return torch.bmm(x, x.permute(0, 2, 1)).view(
            x.shape[0],
            x.shape[1],
            x.shape[1])

    def calculate_coordinates(self, matrix) -> torch.Tensor:
        return self.rnaho.calculate_coordinates_from_distance_matrix(matrix)

    def plot_predicted_structure(self, matrix) -> None:
        self.rnaho.plot_reconstructed_from_dot_molecule(matrix) 


# model = RNAStructureTransformer(number_of_layers=16, d_model=256)
model = ArmnetLikeModel(depth=8, dim=192).to(struc_device) #dim=192
number_of_epochs = 10 # 5
criterion = MaskedLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-5)
# scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer,
#                                                        T_max=number_of_epochs)
# scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.99)
scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=10000, eta_min=1e-6)
# scheduler = torch.optim.lr_scheduler.OneCycleLR(optimizer,
#                                                 max_lr=1e-3,
#                                                 steps_per_epoch=len(train_loader),
#                                                 epochs=number_of_epochs)

# device = 'cuda:0' if torch.cuda.is_available() else 'cpu'


history_train_loss = []
history_val_loss = []

model, hist_train, hist_val = run_experiment(
    model, 
    train_loader,
    val_loader,
    number_of_epochs,
    criterion,
    optimizer,
    scheduler,
    history_train_loss,
    history_val_loss,
)


with open('final_results.pkl', 'wb') as file:
    pickle.dump([model, history_train_loss, history_val_loss], file=file)
