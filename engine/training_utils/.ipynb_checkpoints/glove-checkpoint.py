from armnet_insides.config import ArmNetConfig
from armnet_insides.model import ArmNet
from armnet_insides.dataset import RNA_Dataset_Test
from armnet_insides.training_utils import DeviceDataLoader

from pathlib import Path

import torch
from torch.utils.data import DataLoader
import pandas as pd


def get_armnet(weights_path: str | Path,
               device: str | torch.device) -> ArmNet:
    model = ArmNet(
        depth=ArmNetConfig.num_encoder_layers,
        num_convs=ArmNetConfig.num_conv_layers,
        adj_ks=ArmNetConfig.conv_2d_kernel_size,
        attn_kernel_size=ArmNetConfig.conv_1d_kernel_size,
        dropout=ArmNetConfig.dropout,
        conv_use_drop1d=ArmNetConfig.conv_1d_use_dropout,
        use_bppm=ArmNetConfig.use_bppm,
    )

    model.load_state_dict(
        torch.load(weights_path, map_location=device)
    )
    model.to(device)
    return model

def preprocess_with_shape_armnet(model: ArmNet,
                                 sequences: list[str],
                                 device: str | torch.device
                     ) -> list[torch.Tensor]:
    
    # max_length = max(len(seq) for seq in sequences)

    df = pd.DataFrame({'sequence_id': [0] * len(sequences),
                       'sequence': sequences,
                       'id_min': [0] * len(sequences),
                       'id_max': [len(seq) for seq in sequences]})

    test_dataset = RNA_Dataset_Test(
        df,
        use_bppm=ArmNetConfig.use_bppm,
        bppm_path=None)

    test_dataloader = DeviceDataLoader(
        DataLoader(
            dataset=test_dataset,
            batch_size=len(sequences),
            num_workers=ArmNetConfig.num_workers,
            persistent_workers=True,
            shuffle=False,
            drop_last=False),
        device=device)

    # armnet_embs = []
    with torch.no_grad():
        x, _ = next(iter(test_dataloader))
        armnet_embs = model(x).cpu()[:, 1:-1]

    # seq_batch = torch.cat(armnet_embs, dim=0)
    # return seq_batch
    return armnet_embs
