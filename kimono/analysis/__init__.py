"""Classes for analysis of structural motifs."""


import os
import networkx as nx
import click as ck

import pandas as pd 

from pathlib import Path

from kimono import SEQUENCE_SAVE_DIR, STRUCTURE_SAVE_DIR, RESULTS_SAVE_DIR
from kimono.motif import StructuralMotif

from pydantic import BaseModel

class MotifAnalysisConfig(BaseModel):
    """Configuration for motif analysis."""
    
    """Data directory."""
    data_dir: Path = None # If not none, this overrides the sequence_path and structure_path 

    """Path to save sequences to."""
    sequence_path: Path = SEQUENCE_SAVE_DIR

    """Path to save structures to."""
    structure_path: Path = STRUCTURE_SAVE_DIR
    structure_format: str = "pdb" 
    structure_database: str = "AF" # AlphaFold, PDB, etc.

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
        config: MotifAnalysisConfig = MotifAnalysisConfig(),
    ) -> None:

        self.dataset_path = config.dataset_path

        self.alphafold_structure_dir = config.alphafold_structure_dir

        self.use_dataset = config.use_dataset

        self.af_model_params = config.af_model_params

        self.af_ignore_fragments    = config.af_ignore_fragments
        self.af_model_version       = config.af_model_version
        self.af_file_extension      = config.af_file_extension
 
      
        

        if self.use_dataset == "dbptm":
            self.dbPTM = self._load_dbptm()


        
        
    
    
    def _load_dbptm(self) -> pd.DataFrame:
        """Load the DBPTM database."""
        return pd.read_csv(
            
        )

    def _load_alphafold(self) -> None:
        """Load an alphafold structure."""

        motifs = {}

        

        for site in self.sites: 
            
            
            pdb_path = self.alphafold_structure_dir / self._get_af_filename(
                acc_id,
            )

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
        


