"""Load data from PhosphositePlus (PSP) dataset."""

import pandas as pd

def load_psp_list(
    path_to_psp: str,
):
    df = pd.read_csv(path_to_psp, sep='\t', header=0)