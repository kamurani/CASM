import logging as log
import os
from pathlib import Path
from typing import Callable, Dict, Generator, List, Optional, Tuple
from urllib.error import HTTPError

import networkx as nx
from tqdm import tqdm

from graphein.ml.conversion import GraphFormatConvertor
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graphs_mp
from graphein.protein.utils import (
    download_alphafold_structure,
    download_pdb,
    download_pdb_multiprocessing,
)
from graphein.utils.utils import import_message

try:
    import torch
    from torch_geometric.data import Data, Dataset, InMemoryDataset
except ImportError:
    import_message(
        "graphein.ml.datasets.torch_geometric_dataset",
        "torch_geometric",
        conda_channel="pyg",
        pip_install=True,
    )

import pandas as pd

from CASM.kin_sub_pairs import get_atp_site_dict, KIN_ATP_COORDS_RMSD

"""
Pairs two graphs together in a single ``Data`` instance.
"""
def pair_data(a: Data, b: Data) -> Data:

    out = Data()
    out.a = a 
    out.b = b 
    return out

"""
TODO: handle errors (e.g. missing values for threshold checking in the dataframe)
within this dataset init, instead of filtering ourselves before we pass the dataframe
to this dataset class.
"""
class KinaseSubstrateDataset(Dataset):
    def __init__(
        self,
        root,
        name, 
        pdb_codes: Optional[List[str]] = None,
        uniprot_ids: Optional[List[str]] = None,

        df: pd.DataFrame = None, 

        kinase_metadata_dict: dict = get_atp_site_dict(
            KIN_ATP_COORDS_RMSD,
        ),
        


        graph_labels: Optional[List[torch.Tensor]] = None,
        node_labels: Optional[List[torch.Tensor]] = None,
        chain_selections: Optional[List[str]] = None,

        graphein_config: ProteinGraphConfig = ProteinGraphConfig(),
        graph_format_convertor: GraphFormatConvertor = GraphFormatConvertor(
            src_format="nx", dst_format="pyg"
        ),
        graph_transformation_funcs: Optional[List[Callable]] = None,
        pdb_transform: Optional[List[Callable]] = None,
        transform: Optional[Callable] = None,
        pre_transform: Optional[Callable] = None,
        pre_filter: Optional[Callable] = None,
        num_cores: int = 16,
        af_version: int = 2,
    ):
        """Dataset class for protein graphs.

        Dataset base class for creating graph datasets.
        See `here <https://pytorch-geometric.readthedocs.io/en/latest/notes/
        create_dataset.html>`__ for the accompanying tutorial.

        :param root: Root directory where the dataset should be saved.
        :type root: str
        :param pdb_codes: List of PDB codes to download and parse from the PDB.
            Defaults to ``None``.
        :type pdb_codes: Optional[List[str]], optional
        :param uniprot_ids: List of Uniprot IDs to download and parse from
            Alphafold Database. Defaults to ``None``.
        :type uniprot_ids: Optional[List[str]], optional
        :param graph_label_map: Dictionary mapping PDB/Uniprot IDs to
            graph-level labels. Defaults to ``None``.
        :type graph_label_map: Optional[Dict[str, Tensor]], optional
        :param node_label_map: Dictionary mapping PDB/Uniprot IDs to node-level
            labels. Defaults to ``None``.
        :type node_label_map: Optional[Dict[str, torch.Tensor]], optional
        :param chain_selection_map: Dictionary mapping, defaults to ``None``.
        :type chain_selection_map: Optional[Dict[str, List[str]]], optional
        :param graphein_config: Protein graph construction config, defaults to
            ``ProteinGraphConfig()``.
        :type graphein_config: ProteinGraphConfig, optional
        :param graph_format_convertor: Conversion handler for graphs, defaults
            to ``GraphFormatConvertor(src_format="nx", dst_format="pyg")``.
        :type graph_format_convertor: GraphFormatConvertor, optional
        :param graph_transformation_funcs: List of functions that consume a
            ``nx.Graph`` and return a ``nx.Graph``. Applied to graphs after
            construction but before conversion to pyg. Defaults to ``None``.
        :type graph_transformation_funcs: Optional[List[Callable]], optional
        :param pdb_transform: List of functions that consume a list of paths to
            the downloaded structures. This provides an entry point to apply
            pre-processing from bioinformatics tools of your choosing. Defaults
            to ``None``.
        :type pdb_transform: Optional[List[Callable]], optional
        :param transform: A function/transform that takes in a
            ``torch_geometric.data.Data`` object and returns a transformed
            version. The data object will be transformed before every access.
            Defaults to ``None``.
        :type transform: Optional[Callable], optional
        :param pre_transform:  A function/transform that takes in an
            ``torch_geometric.data.Data`` object and returns a transformed
            version. The data object will be transformed before being saved to
            disk. Defaults to ``None``.
        :type pre_transform: Optional[Callable], optional
        :param pre_filter:  A function that takes in a
            ``torch_geometric.data.Data`` object and returns a boolean value,
            indicating whether the data object should be included in the final
            dataset. Optional, defaults to ``None``.
        :type pre_filter: Optional[Callable], optional
        :param num_cores: Number of cores to use for multiprocessing of graph
            construction, defaults to ``16``.
        :type num_cores: int, optional
        :param af_version: Version of AlphaFoldDB structures to use,
            defaults to ``2``.
        :type af_version: int, optional
        """


        """""""""""""""""""""""""""""""""""""""""" 
        self.name = name



        
        """
        Create ``data_list`` using dataframe 
        """

        # TODO: download all necessary kinases
        # get list of all kinases (can't use unique as KIN_ACC_ID will have list of kinases)

        self.substrates = list(df.SUB_ACC_ID.unique())
        self.kinases = list(df.KIN_ACC_ID.explode().unique())
        
        self.kinase_metadata_dict = kinase_metadata_dict


        # Generate pairs (structures) with labels
        
        

        # have to download pdb files for self.structures; but have to create graphs separately for sub / kin
        self.structures = self.substrates + self.kinases



        self.pairs = [] 

        # TODO generate pairs 

        df2: dict = dict([((sub, mod_rsd), (pos, neg)) for sub, mod_rsd, pos, neg in zip(df.SUB_ACC_ID, df.SUB_MOD_RSD, df.KIN_ACC_ID, df.NEG_KIN_ACC_ID) ])

        # For each site
        for (sub, mod_rsd), (pos, neg) in df2.items():
            
            # Iterate through positive examples
            label = 1
            for kin in pos:
                pair = {
                    "sub": sub,
                    "mod_rsd": mod_rsd,
                    "kin": kin,
                    "label": label,     
                }
                self.pairs.append(pair)    
            # Iterate through negative examples    
            label = 0 
            for kin in neg:     
                pair = {
                    "sub": sub,
                    "mod_rsd": mod_rsd,
                    "kin": kin,
                    "label": label,     
                }
                self.pairs.append(pair)
            
        print(self.pairs[0:10])

        
        self.data_list: List[Tuple[Data, Data]] = None

        """
        
        """


        self.pdb_codes = (
            [pdb.lower() for pdb in pdb_codes]
            if pdb_codes is not None
            else None
        )
        self.uniprot_ids = (
            [up.upper() for up in uniprot_ids]
            if uniprot_ids is not None
            else None
        )

        if self.pdb_codes and self.uniprot_ids:
            self.structures = self.pdb_codes + self.uniprot_ids
        elif self.pdb_codes:
            self.structures = pdb_codes
        elif self.uniprot_ids:
            self.structures = uniprot_ids
        self.af_version = af_version

        # Labels & Chains

        self.examples: Dict[int, str] = dict(enumerate(self.pairs))

        if graph_labels is not None:


            self.graph_label_map = dict(enumerate(graph_labels))

        else:
            self.graph_label_map = None

        if node_labels is not None:
            self.node_label_map = dict(enumerate(node_labels))
        else:
            self.node_label_map = None
        
        if chain_selections is not None:
            self.chain_selection_map = dict(enumerate(chain_selections))
        else:
            self.chain_selection_map = None
        self.validate_input()
        self.bad_pdbs: List[str] = []
        self.bad_uniprot_ids: List[str] = []
        # Configs
        self.config = graphein_config
        self.graph_format_convertor = graph_format_convertor
        self.num_cores = num_cores
        self.pdb_transform = pdb_transform
        self.graph_transformation_funcs = graph_transformation_funcs
        super().__init__(
            root,
            transform=transform,
            pre_transform=pre_transform,
            pre_filter=pre_filter,
        )
        self.config.pdb_dir = Path(self.raw_dir)

    @property
    def raw_file_names(self) -> List[str]:
        """Names of raw files in the dataset."""
        return [f"{pdb}.pdb" for pdb in self.structures]

    @property
    def processed_file_names(self) -> List[str]:
        """Names of processed files to look for"""
        if self.chain_selection_map is not None:
            return [
                f"{pdb}_{chain}.pt"
                for pdb, chain in zip(
                    self.structures, self.chain_selection_map.values()
                )
            ]
        else:
            return [f"{pdb}.pt" for pdb in self.structures]

    def validate_input(self):
        if self.graph_label_map is not None:
            assert len(self.structures) == len(
                self.graph_label_map
            ), "Number of proteins and graph labels must match"
        if self.node_label_map is not None:
            assert len(self.structures) == len(
                self.node_label_map
            ), "Number of proteins and node labels must match"

        # Chain selection
        if True:
            assert len(self.structures) == len(
                self.chain_selection_map
            ), "Number of proteins and chain selections must match"
            assert len(
                {
                    f"{pdb}_{chain}"
                    for pdb, chain in zip(
                        self.structures, self.chain_selection_map
                    )
                }
            ) == len(self.structures), "Duplicate protein/chain combinations"

    def download(self):
        """Download the PDB files from RCSB or Alphafold."""
        self.config.pdb_dir = Path(self.raw_dir)
        if self.pdb_codes:
            # Only download undownloaded PDBs
            to_download = [
                pdb
                for pdb in set(self.pdb_codes)
                if not os.path.exists(Path(self.raw_dir) / f"{pdb}.pdb")
            ]
            
            download_pdb_multiprocessing(
                to_download,
                self.raw_dir,
                max_workers=self.num_cores,
                strict=False,
            )
            self.bad_pdbs = self.bad_pdbs + [
                pdb
                for pdb in set(self.pdb_codes)
                if not os.path.exists(Path(self.raw_dir) / f"{pdb}.pdb")
            ]
        if self.uniprot_ids:

            # Only download undownloaded Uniprot IDs
            to_download = [
                uniprot
                for uniprot in set(self.uniprot_ids)
                if not os.path.exists(Path(self.raw_dir) / f"{uniprot}.pdb")
            ]
            print(f"Uniprot IDs to download: {to_download}")

            for uniprot in tqdm(to_download):

                
                fn = download_alphafold_structure(
                    uniprot,
                    out_dir=self.raw_dir,
                    version=self.af_version,
                    aligned_score=False,
                    rename=True,
                )
                
            self.bad_uniprot_ids = self.bad_uniprot_ids + [
                uniprot
                for uniprot in set(self.uniprot_ids)
                if not os.path.exists(Path(self.raw_dir) / f"{uniprot}.pdb")
            ]
                
                    
            print(f"Bad uniprot ids: {self.bad_uniprot_ids}")    
        
        # TODO: remove bad uniprot / pdb ids from self.structures

    def len(self) -> int:
        """Returns length of data set (number of structures)."""
        return len(self.structures)

    def transform_pdbs(self):
        """
        Performs pre-processing of PDB structures before constructing graphs.
        """
        structure_files = [
            f"{self.raw_dir}/{pdb}.pdb" for pdb in self.structures
        ]
        for func in self.pdb_transform:
            func(structure_files)

    def transform_graphein_graphs(self, graph: nx.Graph):
        for func in self.graph_transformation_funcs:
            graph = func(graph)
        return graph

    def process(self):
        """Processes structures from files into PyTorch Geometric Data."""       
        # Preprocess PDB files
        if self.pdb_transform:
            self.transform_pdbs()


        # Generate graphs for substrates
        #TODO

        torch.load()

        # Generate graphs for kinases 
        #TODO

        d = self.kinase_metadata_dict


        
        coords = d[k]

        idx = 0
        # Chunk dataset for parallel processing
        chunk_size = 128

        def divide_chunks(l: List[str], n: int = 2) -> Generator:
            for i in range(0, len(l), n):
                yield l[i : i + n]

        # chunks = list(divide_chunks(self.structures, chunk_size))
        chunks: List[int] = list(
            divide_chunks(list(self.examples.keys()), chunk_size)
        )

        for chunk in tqdm(chunks):
            pdbs = [self.examples[idx] for idx in chunk]
            # Get chain selections
            if self.chain_selection_map is not None:
                chain_selections = [
                    self.chain_selection_map[idx] for idx in chunk
                ]
            else:
                chain_selections = ["all"] * len(chunk)

            # Create graph objects
            file_names = [f"{self.raw_dir}/{pdb}.pdb" for pdb in pdbs]
            graphs = construct_graphs_mp(
                pdb_path_it=file_names,
                config=self.config,
                chain_selections=chain_selections,
                return_dict=False,
            )
            if self.graph_transformation_funcs is not None:
                graphs = [self.transform_graphein_graphs(g) for g in graphs]

            # Convert to PyTorch Geometric Data
            graphs = [self.graph_format_convertor(g) for g in graphs]

            # Assign labels
            if self.graph_label_map:
                labels = [self.graph_label_map[idx] for idx in chunk]
                for i, _ in enumerate(chunk):
                    graphs[i].graph_y = labels[i]
            if self.node_label_map:
                labels = [self.node_label_map[idx] for idx in chunk]
                for i, _ in enumerate(chunk):
                    graphs[i].graph_y = labels[i]

            data_list = graphs

            del graphs

            if self.pre_filter is not None:
                data_list = [g for g in data_list if self.pre_filter(g)]

            if self.pre_transform is not None:
                data_list = [self.pre_transform(data) for data in data_list]

            for i, (pdb, chain) in enumerate(zip(pdbs, chain_selections)):
                if self.chain_selection_map is None:
                    torch.save(
                        data_list[i],
                        os.path.join(self.processed_dir, f"{pdb}.pt"),
                    )
                else:
                    torch.save(
                        data_list[i],
                        os.path.join(self.processed_dir, f"{pdb}_{chain}.pt"),
                    )

    def get(self, idx: int):
        """
        Returns PyTorch Geometric Data object for a given index.

        :param idx: Index to retrieve.
        :type idx: int
        :return: PyTorch Geometric Data object.
        """
        # TODO: work out structure for storing the pairs and labels
        
        # structures[idx]['kinase'] # kin, sub, mod_rsd, label, 
        # OR alternative structure: 

        # Load kinase, substrate, mod_rsd, label

        kinase = torch.load(
                os.path.join(self.processed_dir, f"KIN_{self.examples[idx]['kin']}.pt")
        )

        site = torch.load(
                os.path.join(self.processed_dir, f"SUB_{self.examples[idx]['sub']}_{self.examples[idx]['mod_rsd']}.pt")
        )

        label = torch.tensor(self.examples[idx]['label'])


        if self.chain_selection_map is not None:
            return torch.load(
                os.path.join(
                    self.processed_dir,
                    f"{self.structures[idx]}_{self.chain_selection_map[idx]}.pt",
                )
            )
        else:

            return torch.load(
                os.path.join(self.processed_dir, f"{self.structures[idx]}.pt")
            )

