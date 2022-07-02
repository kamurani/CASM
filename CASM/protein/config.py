"""Base config object for use with protein clustering."""
# CASM
# Author: C. I. McM <c.mcmenamie@unsw.edu.au>
# License: MIT
# Project Website: https://github.com/cimranm/CASM
# Code Repository: https://github.com/cimranm/CASM

from typing import List, Union, Callable, Optional
from pydantic import BaseModel


class ClusterConfig(BaseModel):

    unsupervised: bool = True
    node_embedding: Optional[List[Union[Callable, str]]] = "meiler"
