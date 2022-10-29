"""Load in phylogenetic tree to calculate distance matrix"""



import csv
from typing import Dict
import dendropy 

import click as c


"""
TODO: make more generalisable by generating a function that returns a dict
for any two column names in the file.  

i.e. similar to graphein "meta programming" thing in that class that generated
functions without pre-coding all of them 

e.g. can say "_xName2Group" and it will generate a function that creates a dict for that file 
as expected.


TODO: dump this into a file in python format, so we no longer have to read in from a file everytime 
we can just import it as a module.
"""
def get_kinase_name_map(
    filepath: str, 
    from_name: str, # Which column name to use as the map input

) -> Dict[str, Dict[str, str]]:

    d = {}
    with open(filepath) as f:
        tsv_file = csv.DictReader(f, delimiter="\t")

        for line in tsv_file:

            new_line = {
                key.strip(): value.strip() for key, value in line.items()
            }
            
            idx = new_line[from_name]
            #print(idx)
            d[idx] = new_line
    
    return d

def get_phylo_tree(
    filepath: str, 
    
):

    tree = dendropy.Tree.get(
            path=filepath, 
            schema="newick",
    )
    return tree

"""
Returns a dictionary used to look up pairwise distances
of taxa (in Uniprot ID format).  Self pairs are 0. 
"""
def get_distance_matrix(
    filepath: str,                  # Path to the tree file
    kinase_table: str,              # Path to the lookup table
) -> Dict[str, Dict[str, float]]:

    tree = get_phylo_tree(filepath)
    pdm = tree.phylogenetic_distance_matrix()

    # Convert to dictionary 
   
    # Go through distinct pairs, and create distance matrix
    matrix = {}

    # Set empty dicts for the first index 
    labels = tree.taxon_namespace.labels()

    from_name = "Manning Name"
    name_map: dict = get_kinase_name_map(filepath=kinase_table, from_name=from_name)

    uniprot_names = [name_map[label.split()[-1]]["UniprotID"] for label in labels] # all last name labels in the tree, as Uniprot IDs

    for a in uniprot_names:
        # Initialise
        matrix[a] = {}

        # Self-distance is zero
        matrix[a][a] = 0


    count = 0
    for (t1, t2) in pdm.distinct_taxon_pair_iter():

        count += 1

        # Get distance val
        distance = pdm.distance(t1, t2)

        # Get last word in label, which will be kinase's Manning Name
        x, y = t1.label.split()[-1], t2.label.split()[-1] 

        # Convert to Uniprot acc ID
        x_uniprot = name_map[x]["UniprotID"]
        y_uniprot = name_map[y]["UniprotID"]

        # Save 2nd dict that stores the distance
        try:
            matrix[x_uniprot][y_uniprot] = distance
            # Save other way round 
            matrix[y_uniprot][x_uniprot] = distance
        except:
            raise Exception(f"Something went wrong with initialising first matrix")

    return matrix

@c.command()
@c.argument(
    'TREE_PATH', nargs=1,
    type=c.Path(exists=True, file_okay=True, dir_okay=False),
)
def main(
    tree_path,
):

    # Tree path in PHYLIP format
    pdc = get_distance_matrix(tree_path)

    t1 = dendropy.Taxon("CAMK MAPKAPK MNK MNK2")
    t2 = dendropy.Taxon("CMGC CDK CRK7 CRK7")
    dists = pdc(t1, t2)
    print(dists)


if __name__ == "__main__":
    main()
