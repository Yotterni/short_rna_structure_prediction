from rna_holder import RNAHolder
from pathlib import Path

import os
import pickle


def prepare_rf_dataset(rf_pdbs_path: str | Path,
                       dataset_path: str | Path) -> None:
    counter = 0
    pdb_filenames = sorted(os.listdir(rf_pdbs_path))
    for idx, filename in enumerate(pdb_filenames):
        rnaho = RNAHolder(os.path.join(rf_pdbs_path, filename))
        onehot, pairwise_dist_matrix = rnaho.get_representation_and_matrix()

        if idx % 10000 == 0:
            print(f'Processed {idx} molecules')

        with open(
                os.path.join(dataset_path,
                             f'{filename.split('.')[0]}__{counter}.pkl'), 'wb') as file:
            pickle.dump((rnaho.sequence, onehot, pairwise_dist_matrix), file)

if __name__ == '__main__':
    prepare_rf_dataset('data/rf_diffusion_data',
                       'data/rf_dataset')
