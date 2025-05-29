from graphein.rna.graphs import construct_graph
from graphein.rna.visualisation import plotly_rna_structure_graph

from networkx.classes.graph import Graph
from pathlib import Path

import typing as tp

import numpy as np
import pandas as pd
import plotly.express as px
import seaborn as sns
import torch
import torch.nn.functional as F


class RNAHolder:
    """
    A universal class for processing RNA molecules in the Stanford
    RNA3D competition.
    """
    def __init__(self, path_to_pdb: tp.Optional[Path | str] = None,
                 cache_graph: bool = False) -> None:
        self.cache_graph = cache_graph
        self.c1s_atoms = []
        self.coordinates = None
        self.sequence = ''

        if path_to_pdb is not None:
            self.graph = construct_graph(path=path_to_pdb)
            self.parse_graph_(self.graph.nodes(data=True))

        if not cache_graph:
            self.graph = None

        self.pairwise_euclidean_distance_matrix = None
        self.pairwise_dot_matrix = None

    def parse_graph_(self, nodes) -> list[dict[str, tp.Any]]:
        """
        From networkx graph to list of c1' atoms.
        :param nodes:
        :return:
        """
        nucleotides = []
        coordinates = []
        for atom_name, atom_features in nodes:
            atom_type = atom_features['atom_type']
            if atom_type == "C1'":
                atom_features.update({'name': atom_name})
                self.c1s_atoms.append(atom_features)
                coordinates.append(
                    torch.Tensor(atom_features['coords']))
                nucleotides.append(atom_features['residue_name'])
        self.sequence = ''.join(nucleotides)
        self.coordinates = torch.stack(coordinates, dim=0)
        return self.c1s_atoms

    def plot_original_pdb_structure(self) -> None:
        """
        Takes the original pdb structure and plots a plotly 3d image.
        :return:
        """
        if not self.cache_graph:
            print("Original structure isn't cached")
            return
        fig = plotly_rna_structure_graph(
            self.graph,
            node_alpha=1,
            node_size_multiplier=1,
            node_size_min=3,
            colour_edges_by="kind",
            colour_nodes_by="seq_position",
            label_node_ids=False,
        )
        fig.show()

    def calculate_pairwise_euclidean_distance_matrix(self,
                                           coordinates: tp.Optional[torch.Tensor] = None
                                           ) -> torch.Tensor:
        """
        Calculates pairwise euclidean_distance matrix.
        :return:
        """
        if coordinates is None:
            coordinates = self.coordinates
        self.pairwise_euclidean_distance_matrix = torch.cdist(
            coordinates, coordinates, p=2).view(
            len(self.coordinates), len(self.coordinates)) # ** 2
        return self.pairwise_euclidean_distance_matrix

    def calculate_pairwise_dot_matrix(self,
            coordinates: tp.Optional[torch.Tensor] = None) -> torch.Tensor:
        if coordinates is None:
            coordinates = self.coordinates

        if len(coordinates.shape) == 2:
            centroid = coordinates.mean(axis=0)[None]
            coordinates -= centroid
            self.pairwise_dot_matrix = coordinates @ coordinates.T
        else:
            centroid = coordinates.mean(axis=1)[:, None, :]
            coordinates -= centroid
            self.pairwise_dot_matrix = torch.bmm(coordinates, coordinates.T)
        return self.pairwise_dot_matrix

    def plot_pairwise_euclidean_distance_matrix(self) -> None:
        """
        Calculates pairwise euclidean_distance matrix.
        :return:
        """
        if self.pairwise_euclidean_distance_matrix is None:
            self.pairwise_euclidean_distance_matrix = (
                self.calculate_pairwise_euclidean_distance_matrix())

        matrix = self.pairwise_euclidean_distance_matrix.detach().cpu().numpy()
        sns.heatmap(matrix, cmap=sns.color_palette(
            "ch:start=.2,rot=-.3", as_cmap=True))

    def plot_pairwise_dot_matrix(self) -> None:
        if self.pairwise_dot_matrix is None:
            self.pairwise_dot_matrix = self.calculate_pairwise_dot_matrix()

        matrix = self.pairwise_dot_matrix.detach().cpu().numpy()
        sns.heatmap(matrix, cmap=sns.color_palette(
            "ch:start=.2,rot=-.3", as_cmap=True))

    def calculate_coordinates_from_euclidean_distance_matrix(
            self,
            matrix: tp.Optional[torch.Tensor] = None
    ) -> np.ndarray:
        """
        Converts pairwise euclidean_distance matrix to coordinates.
        :return:
        """
        if matrix is None:
            dij_matrix = self.pairwise_euclidean_distance_matrix ** 2
        else:
            dij_matrix = matrix ** 2
        d1j_matrix = dij_matrix[0].view(1, -1)
        d1j_matrix = d1j_matrix.repeat(dij_matrix.shape[0], 1)

        di1_matrix = dij_matrix[:, 0].view(-1, 1)
        di1_matrix = di1_matrix.repeat(1, dij_matrix.shape[1])

        mij_matrix = (1 / 2) * (d1j_matrix + di1_matrix - dij_matrix)

        u, s, v = np.linalg.svd(mij_matrix, full_matrices=False)
        u = u[:, :3]
        return u * np.sqrt(s[:3]).reshape(-1, 1).T

    def plot_reconstr_mol_from_eucl_distance_matrix(
        self,
        euclidean_distance_matrix: tp.Optional[torch.Tensor] = None) -> None:
        """
        Plots reconstructed molecule.
        :return:
        """
        if euclidean_distance_matrix is None:
            euclidean_distance_matrix = self.pairwise_euclidean_distance_matrix

        coordinates = (
            self.calculate_coordinates_from_euclidean_distance_matrix(euclidean_distance_matrix))
        df = pd.DataFrame(coordinates)
        df.columns = ['x', 'y', 'z']
        img = px.scatter_3d(df, x='x', y='y', z='z',
                            title='reconstructed molecule',
                            color=df.index)
        img.show()

    def calculate_coordinates_from_dot_distance_matrix(
            self,
            dot_matrix: tp.Optional[torch.Tensor] = None) -> torch.Tensor:
        if dot_matrix is None:
            dot_matrix = self.pairwise_dot_matrix

        dot_matrix = dot_matrix.cpu()

        u, s, _ = np.linalg.svd(dot_matrix, full_matrices=False)
        u = u[:, :3]
        return u * np.sqrt(s[:3]).reshape(-1, 1).T
    
    def plot_reconstructed_from_dot_molecule(self, dot_matrix = None) -> None:
        coordinates = (
            self.calculate_coordinates_from_euclidean_distance_matrix(dot_matrix))
        # print(coordinates.shape)
        df = pd.DataFrame(coordinates)
        df.columns = ['x', 'y', 'z']
        # print(df.shape)
        img = px.scatter_3d(df, x='x', y='y', z='z',
                   title='reconstructed molecule',
                   color=df.index)
        img.show()
    
    def calculate_onehot_embedding(self, sequence: tp.Optional[str] = None
        ) -> torch.Tensor:
        """
        Calculates onehot embedding for the sequence.
        :return:
        """
        if sequence is None:
            sequence = self.sequence
        nt_to_ind = {'A': 0, 'U': 1, 'G': 2, 'C': 3}
        nt_indexes = torch.tensor([nt_to_ind[nt] for nt in sequence])
        return F.one_hot(nt_indexes, num_classes=4).to(torch.float)

    def get_representation_and_dot_matrix(self
                                      ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Just calls `self.calculate_onehot_embedding()` and
        `self.calculate_pairwise_dot_matrix()`.
        :return:
        """
        return (self.calculate_onehot_embedding(),
                self.calculate_pairwise_dot_matrix())


    def get_representation_and_eucl_matrix(self
                                      ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Just calls `self.calculate_onehot_embedding()` and
        `self.calculate_pairwise_euclidean_distance_matrix()`.
        :return:
        """
        return (self.calculate_onehot_embedding(),
                self.calculate_pairwise_euclidean_distance_matrix())
