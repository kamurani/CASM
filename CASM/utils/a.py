



"Get AlphaFold filename from database download"
def get_pdb_filename(
    acc_id: str, 
    file_extension: str, # could be .pdb.gz 

    model_version: int = 3, 
    ignore_fragments: bool = True, 

) -> str: 

    if ignore_fragments:
        fragment = 1
        return f"AF-{acc_id}-F{fragment}-model_v{model_version}.{file_extension}"
    else:
        raise NotImplementedError(f"Multiple fragments for AF structure not implemented. ")

def is_list_of_str(
    lst,
):
    if lst and isinstance(lst, list):
        return all(isinstance(elem, str) for elem in lst)
    else: 
        return False
