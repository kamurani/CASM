"""Get binding site locations for kinases"""

import os
from pathlib import Path
from typing import Callable, List, Union

import click as c

import json

import requests


class AlphaFillLoader():

    def __init__(

        self,
        
        cif_dir: Union[str, Path],
        structure_extension: str = "cif", 
        metadata_extension: str = "json",
        config = None, # config object for graph construction 
        ignore: bool = False, # ignore missing directories 

        functions: List[Callable] = None, # TODO: specify subgraph selecting functions to be applied 
    ):



        # TODO: store `pdb_dir` in the config , more elegant
        # TODO: check if exists, and contains expected filename extensions
        if type(cif_dir) is str: 
            cif_dir = Path(cif_dir)

        self.cif_dir = cif_dir
        if not ignore and not os.path.isdir(str(cif_dir)):
            raise ValueError(f"No such directory '{str(cif_dir)}'") 

        self.config = config
        self.structure_extension = structure_extension
        self.metadata_extension = metadata_extension

    def _structure_url(self, acc_id: str):
        
        return f"https://alphafill.eu/v1/aff/{acc_id}-F1"

    def _metadata_url(self, acc_id: str):

        return f"https://alphafill.eu/v1/aff/{acc_id}-F1/json"

    def _get_filename(self, acc_id: str, file_extension: str):
        return f"{acc_id}-ALPHAFILL.{file_extension}"
    

    def _structure_path(self, acc_id: str):
        
        p = self.cif_dir / self._get_filename(acc_id, self.structure_extension)
        p = p.resolve() # get absolute path
        path = str(p)
        return path

    def _metadata_path(self, acc_id: str):
        
        p = self.cif_dir / self._get_filename(acc_id, self.metadata_extension)
        p = p.resolve() # get absolute path
        path = str(p)
        return path




    """
    Download 
    """
    def download_alphafill(
        self,
        acc_id: str,
        overwrite: bool = True, # TODO # Overwrite already existing files with same name
    ):
        try:
            # Download metadata
            response = requests.get(self._metadata_url(acc_id))
            open(self._metadata_path(acc_id), "wb").write(response.content)

            # Download structure
            response = requests.get(self._structure_url(acc_id))
            open(self._structure_path(acc_id), "wb").write(response.content)
        except: 
            print(f"Failed to download {acc_id} from AlphaFill.")

    def download_alphafill_list(
        self,
        acc_id_list: List[str],
    ):

        for acc in acc_id_list:
            self.download_alphafill(acc)

    """
    Parses AlphaFill metadata JSON file to find the best (lowest RMSD) scoring 
    transplant of a specified molecule ``label_comp_id``.

    Uses ``local 
    """
    def get_asym_id_lowest_rmsd(
        self,
        acc_id: str, 
        molecule: str = "ATP",

    ):

        # Check existence of path
        fn = self._metadata_path(acc_id)

        if not os.path.isfile(fn): raise FileNotFoundError(f"No file found: {fn}")

        # Parse the JSON file 
        with open(fn) as f:

            json.load(f)



    """
    Get coordinates of a particular atom for a file
    """
    def get_coordinates(
        self,
        acc_id: str,
        
    ):
        # TODO
        path = self._get_structure_path(acc_id)

        # 

        #molecule = 
        (x,y,z) = (0,0,0)
        print(f"{acc_id}\t{x}\t{y}\t{z}")
        

"""
Download 
"""
def get_kinase_binding_site(
    kin_acc_id_list: List[str], 
):

    """
    TODO: convert isoforms to original for AF2 models
    e.g. 

    Q9BUB5-2 --> Q9BUB5
    """


def load_from_file(
    filepath: str = "",
) -> List[str]:

#    Load in from file (list of kinases)
    filepath = "datasets/all_kin_acc_ids.txt"
    # TODO: store default somewhere else


    with open(filepath, 'r', encoding='utf-8') as infile:

        acc_id_list = [line.strip() for line in infile]

    return acc_id_list




@c.command()
@c.argument('ALPHAFILL_DIR', nargs=1,
    type=c.Path(exists=True, file_okay=False, dir_okay=True),
)
def main(
    alphafill_dir,
):

    print(f"Saving JSON and CIF files to {alphafill_dir}")
    acc_ids = load_from_file()

    #acc_ids = ["Q9UHD2", "Q9WTU6"]

    # Alphafill




    loader = AlphaFillLoader(
        cif_dir=alphafill_dir,
    )

    # Download (if required)
    #loader.download_alphafill_list(acc_ids)

    loader.get_coordinates()



if __name__ == "__main__":
    main()

