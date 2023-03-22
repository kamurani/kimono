"""Util functions"""

from kimono.utils.definitions import NODE_ID_STR
from kimono.protein.data import protein_letters_1to3


def get_node_id_string(
    chain_id: str,
    residue: str,
    pos: int,
) -> str:
    """Get a string representation of a node."""
    # TODO: break if not a valid AA
    if len(residue) == 3:
        residue = residue.upper() 
    elif len(residue) == 1:
        residue = protein_letters_1to3[residue.upper()]
    else:
        raise ValueError(f"Invalid residue: {residue}")
    
    return NODE_ID_STR.format(chain_id.upper(), residue, str(pos))