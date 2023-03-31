"""Classes and functions for protein post-translational modifications."""

from kimono.protein.data import protein_letters_3to1
from kimono.utils.utils import get_node_id_string

class PTMSite():
    """A site of a post-translational modification on a protein."""

    symbol_map = {
        "phosphorylation": "p",
        "ubiquitination": "u",
        "acetylation": "ac",
        "methylation": "me",
    }

    def __init__(
        self,
        entry_name: str, 
        acc_id: str,
        residue: str, # 1-letter code 
        position: int,
        ptm_type: str,
        chain_id: str = 'A', 
    ) -> None:

        self.entry_name: str = entry_name # UniProt entry name
        self.acc_id: str = acc_id # uniprot ID that the PTM site is mapped to.

        self.residue: str = residue 
        self.position: int = position
        self.ptm_type: str = ptm_type 
        self.chain_id: str = chain_id

        if self.chain_id !='A':
            raise NotImplementedError("Only chain A is supported for now.")

        if len(self.residue) == 3:
            self.residue = protein_letters_3to1[self.residue.capitalize()] 
        elif len(self.residue) == 1:
            self.residue = self.residue.upper()
        else:
            raise ValueError(f"Invalid residue: {self.residue}")

        self._node_id: str = get_node_id_string(self.chain_id, self.residue.upper(), self.position)

    def __repr__(self):
        return f"PTMSite({self.acc_id}, {self.node_id}-{PTMSite.symbol_map[self.ptm_type.lower()]})" 

    @property
    def node_id(self) -> str:
        self._node_id = get_node_id_string(self.chain_id, self.residue.upper(), self.position)
        return self._node_id 

    


        