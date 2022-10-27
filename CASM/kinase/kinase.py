"""Get binding site locations for kinases"""

import os
from pathlib import Path
from typing import Callable, Dict, List, Union

import click as c

import json

import requests


# CIF file parsing
from gemmi import cif


class AlphaFillLoader():

    def __init__(

        self,
        
        cif_dir: Union[str, Path],
        structure_extension: str = "cif", 
        metadata_extension: str = "json",
        config = None, # config object for graph construction 
        ignore: bool = False, # ignore missing directories 

        functions: List[Callable] = None, # TODO: specify subgraph selecting functions to be applied 
        verbose: bool = True,
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

        self.verbose = verbose

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
        if self.verbose: print(f"Saving JSON and CIF files to {self.cif_dir}")
        for acc in acc_id_list:
            self.download_alphafill(acc)

    """
    Parses AlphaFill metadata JSON file to find the best (lowest RMSD) scoring 
    transplant of a specified molecule ``label_comp_id``.

    Uses local RMSD for backbone within 6Å of the transplanted molecule. 

    Returns dictionary of 
        analogue_id: 
        asym_id: 
        rmsd:
    """
    def get_asym_id_hits(
        self,
        acc_id: str, 
        molecule: str = "ATP",

    ) -> Dict[str, str]:

        # Check existence of path
        fn = self._metadata_path(acc_id)

        if not os.path.isfile(fn): raise FileNotFoundError(f"No file found: {fn}")

        # Parse the JSON file 
        with open(fn) as f:

            data = json.load(f)


            # OK, let's attack this JSON and find the lowest RMSD for our molecule
            
            try: 
                hits: list = data['hits']
            except:
                return None

            # Get all hits for our molecule    
            candidates = []
            for hit in hits:

                for t in hit['transplants']:

                    if t['analogue_id'] == molecule:

                        c = dict(
                            analogue_id=t['analogue_id'],
                            asym_id=t['asym_id'],
                            rmsd=t['rmsd'],
                        )
                        candidates.append(c)

            # Sort by RMSD, increasing order
            c_list = sorted(candidates, key=lambda d: d['rmsd'])

            # Return object with lowest RMSD
            return c_list


    """
    Gets ``asym_id`` of lowest scoring transplant molecule for a list of 
    AlphaFill models 
    """
    def get_asym_id_list(
        self,
        acc_id_list: List[str], 
        molecule: str = "ATP",
    ) -> Dict[str, Dict]:

        out = {}
        for acc_id in acc_id_list:

            try: 
                data = self.get_asym_id_lowest_rmsd(acc_id=acc_id, molecule=molecule)

            except:
                # If no result found, set to empty dictionary
                data = {}

            out[acc_id] = data 

        return out

    """
    Get coordinates of a particular atom for a file

    e.g.

    Gamma phosphorus for ATP (at kinase binding site)
    """
    def get_coordinates(
        self,
        acc_id: Union[str, List[str]],
        

        rmsd_threshold: float = 6.0, # Ignore any data that has RMSD above this value # TODO

        label_asym_id: str = None,      # This needs to be specified if there are multiple transplants of same molecule 
                                        # By default this will use the metadata in same ``cif_dir`` directory from AlphaFill
        type_symbol: str = "P",
        label_atom_id: str = "PG",
        label_comp_id: str = "ATP",
        
    ):
        if type(acc_id) == str:
            # Check existence of path
            fn_cif = self._structure_path(acc_id)
            if not os.path.isfile(fn_cif): raise FileNotFoundError(f"No file found: {fn_cif}")

            # Get ``asym_id`` of lowest scoring molecule (if not already specified to the function)

            if not label_asym_id:

                # try / except here TODO
                try:
                    hits = self.get_asym_id_hits(acc_id=acc_id, molecule=label_comp_id)

                    data: dict = hits[0]

                    label_asym_id   = data['asym_id']
                    rmsd            = data['rmsd']
                except:
                    return dict(coords=None, rmsd=-1)

            # Parse CIF file
            doc = cif.read_file(fn_cif)
            block = doc.sole_block() 

            table = block.find(
                "_atom_site.", 
                ["type_symbol", "label_atom_id", "label_comp_id", "label_asym_id",  # Used in the query checking (input)
                "Cartn_x", "Cartn_y", "Cartn_z"]                                    # cartesian coordinates      (output)
            )
            tags = table.tags 

            query = [type_symbol, label_atom_id, label_comp_id, label_asym_id]
            # TODO more elegant way of checking row match?


            # TODO: get ALL molecules (e.g. all ATPs), and then keep them in a list along with 


            # iterate through rows
            coords = {}
            for row in table:
                
                if (
                    row["_atom_site.type_symbol"]   == type_symbol      and
                    row["_atom_site.label_atom_id"] == label_atom_id    and
                    row["_atom_site.label_comp_id"] == label_comp_id    #and
                    #row["_atom_site.label_asym_id"] == label_asym_id


                    # TODO: take into account missing PG atom (i.e. ATP analogue without phosphate group; still use
                    # a diff atom's coordinates?)
                ):

                    # Save match, indexing using ASYM ID
                    coords[row["_atom_site.label_asym_id"]] = (
                        row["_atom_site.Cartn_x"],
                        row["_atom_site.Cartn_y"],
                        row["_atom_site.Cartn_z"],
                    )

            # no match
            
            if not coords: 
                if self.verbose: print(f"No element found for query: {query}")#raise Exception(f"No element found for query: {query}")
                
        

            
            # Select lowest RMSD. (or lowest that works at least.)
            # i.e. go through hits in ascending rmsd order, until we find one that has coords
            for h in hits:
                asym_id = h['asym_id']
                if asym_id in coords: 
                    return dict(
                        coords=coords[asym_id], 
                        rmsd=rmsd, label_atom_id=label_atom_id,
                        label_comp_id=label_comp_id,
                        label_asym_id=label_asym_id,
                    )


            # Failed
            return dict(coords=None, rmsd=-1)        

       

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

    
    acc_ids = load_from_file()

    #acc_ids = ["Q9UHD2", "Q9WTU6"]

    #acc_ids = ["Q96PY6"] # TODO: for some reason, this alphafill CIF file doesn't have all the atoms for ATP 'Y' ??

    # Alphafill
    loader = AlphaFillLoader(
        cif_dir=alphafill_dir,
    )

    # Download (if required)
    #loader.download_alphafill_list(acc_ids)


    for acc_id in acc_ids:

        try:
            
            data = loader.get_coordinates(
                acc_id=acc_id,
            )
            (x, y, z)   = data['coords']
            rmsd        = data['rmsd']

            molecule    = data["label_comp_id"]
            asym_id     = data["label_asym_id"]
            atom_id     = data["label_atom_id"]
            l = [acc_id, rmsd, atom_id, molecule, asym_id, x, y, z, ]

        #try:
        #    pass  
        except:

            l = [acc_id, 'UNKNOWN']

        print('\t'.join([str(x) for x in l]))

    return


    ATP: dict = loader.get_asym_id_list(acc_ids)

    # print
    for k in ATP.keys():
        try:
            print(f"{k}\t{ATP[k]['asym_id']}\t{ATP[k]['rmsd']}")
        except:
            print(k, "NO DATA")

if __name__ == "__main__":
    main()

