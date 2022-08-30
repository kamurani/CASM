"""Utilities for amino acid sequence"""


from graphein.utils.utils import protein_letters_3to1_all_caps as aa3to1


#aa3to1




def aa1to3(
    aa: str, 
) -> str:
    
    protein_letters_1to3 = {
        "A": "Ala",
        "C": "Cys",
        "D": "Asp",
        "E": "Glu",
        "F": "Phe",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "K": "Lys",
        "L": "Leu",
        "M": "Met",
        "N": "Asn",
        "P": "Pro",
        "Q": "Gln",
        "R": "Arg",
        "S": "Ser",
        "T": "Thr",
        "V": "Val",
        "W": "Trp",
        "Y": "Tyr",
    }

    try: 
        return map[aa.upper()].upper()
    except:
        raise ValueError(f"Specified residue '{aa}' invalid.")