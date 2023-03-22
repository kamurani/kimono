"""Classes for analysis of structural motifs."""


import os
from typing import List
import networkx as nx
import click as ck

import pandas as pd 

from pathlib import Path

from kimono import SEQUENCE_SAVE_DIR, STRUCTURE_SAVE_DIR, RESULTS_SAVE_DIR
from kimono.motif import StructuralMotif

from kimono.protein.data import protein_letters_1to3, protein_letters_3to1

from kimono.utils.utils import get_node_id_string


from pydantic import BaseModel

class PTMSite():
    """A site of a post-translational modification on a protein."""
    def __init__(
        self,
        acc_id: str,
        residue: str, # 1-letter code 
        position: int,
        ptm_type: str,
        chain_id: str = 'A', 
    ) -> None:

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

        self.node_id: str = get_node_id_string(self.chain_id, self.residue, self.position)

    def __repr__(self):
        return f"PTMSite({self.acc_id}, {self.node_id})" 


        





class MotifAnalysisConfig(BaseModel):
    """Configuration for motif analysis."""


    
    """Data directory."""
    data_dir: Path = None # If not none, this overrides the sequence_path and structure_path 

    dbptm_path: Path = data_dir / "dbptm" / "dbptm.csv" if data_dir is not None else None

    """Path to save sequences to."""
    sequence_path: Path = SEQUENCE_SAVE_DIR

    """Path to save structures to."""
    structure_path: Path = STRUCTURE_SAVE_DIR
    structure_format: str = "pdb" 
    structure_database: str = "AF" # AlphaFold, PDB, etc.

    isoform_path: Path = None # Unused for now
    
    dataset_path: Path = None 

    alphafold_structure_dir: Path = None

    af_ignore_fragments: bool = True
    af_model_version: int = 3
    af_file_extension: str = "pdb.gz"

    af_model_params: dict = {
        "ignore_fragments": True,
        "model_version": 3, 
        "file_extension": "pdb.gz", 
    }
    

    """Which PTM dataset to use."""
    use_dataset: str = "dbPTM" # dbptm, ptmnet, etc.

    """Path to save intermediate results to."""
    result_path: Path = RESULTS_SAVE_DIR

    """Settings"""
    force_download: bool = False # If True, will download all data from scratch regardless of whether it already exists.

    def __init__(self, **data):
        super().__init__(**data)

        self.structure_path = Path(self.structure_path) / self.structure_database.lower()

        force = False
        if not self.structure_path.exists():
            if force is True:
                os.makedirs(self.structure_path)
            else:
                raise ValueError(f"Directory {self.structure_path} does not exist.")

        if self.isoform_path is None:
            pass 

    def __repr__(self):
        return f"MotifAnalysisConfig(data_dir={self.data_dir}, sequence_path={self.sequence_path}, structure_path={self.structure_path}, structure_format={self.structure_format}, structure_database={self.structure_database}, result_path={self.result_path}, force_download={self.force_download})"


class MotifAnalysis():

    def __init__(
        self,
        config: MotifAnalysisConfig = None,
    ) -> None:

        self.dataset_path = config.dataset_path

        self.alphafold_structure_dir = config.alphafold_structure_dir

        self.use_dataset = config.use_dataset

        self.af_model_params = config.af_model_params

        self.af_ignore_fragments    = config.af_ignore_fragments
        self.af_model_version       = config.af_model_version
        self.af_file_extension      = config.af_file_extension
 
        self.dbptm_path = config.dbptm_path

        self.sites = []

        if self.use_dataset == "dbptm":
            self.dataset = self._load_dbptm()

        
        
        


        
        
    
    
    def _load_dbptm(self) -> pd.DataFrame:
        """Load the DBPTM database."""
        pass
        self.df = pd.read_csv(
            self.dbptm_path, 
        )

        # Create list of sites for each PTM site recorded in the dataset. 


        site = PTMSite(

        )




    def _load_alphafold(self) -> None:
        """Load an alphafold structure."""

        motifs = {}

        self.sites: List[PTMSite]       # set of sites (i.e. ID and residue position)
        self.motifs                     # subgraph centred around each site # TODO make 
                                        # this efficient without multiple subgraphs; but 
                                        # rather just a refernce to the original graph 

        for site in self.sites: 
            
            pdb_path = self.alphafold_structure_dir / self._get_af_filename(site.acc_id)

            

        return motifs 



    def _get_af_filename(
        self,
        acc_id: str, 
    ) -> str: 

        if self.ignore_fragments:
            fragment = 1
            return f"AF-{acc_id}-F{fragment}-model_v{self.model_version}.{self.file_extension}"
        else:
            raise NotImplementedError(f"Multiple fragments for AF structure not implemented.")
    

    """
    Analysis methods
    """
    def difference_transform(
        self,
    ) -> None:
        """Perform a difference transform on the set of motifs."""
        


