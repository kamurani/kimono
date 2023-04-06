"""Classes for analysis of structural motifs."""


import os
from typing import List
import networkx as nx
import click as ck

from tqdm import tqdm

import pandas as pd 

from pathlib import Path

from kimono import SEQUENCE_SAVE_DIR, STRUCTURE_SAVE_DIR, RESULTS_SAVE_DIR
from kimono.motif import StructuralMotif

from kimono.protein.data import protein_letters_1to3, protein_letters_3to1

from kimono.utils.utils import get_node_id_string

from kimono.ptm import PTMSite


from pydantic import BaseModel







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
    
    max_sites: int = 1000 # Maximum number of PTM sites to include in the analysis.

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

    radius: float = 12.0 # Radius of motif subgraph in Ångströms

    

    """Which PTM dataset to use."""
    use_dataset: str = "dbPTM" # dbptm, ptmnet, etc.

    """Path to save intermediate results to."""
    result_path: Path = RESULTS_SAVE_DIR

    """Settings"""
    force_download: bool = False # If True, will download all data from scratch regardless of whether it already exists.

    include_isoforms: bool = False # If True, will include sequence isoforms in the analysis.

    species_filter: List[str] = None # If not None, will only include sequences from the given species.


    """Filter dict for storing a parameter (i.e. column in dataset) and a threshold value"""
    # plddt score (if available); need to use a file first from a given path in config to then join;
    # to the dataframe.  If no path exists (none) and threshold specified; then the class automatically 
    # loads from AF files as we go. 

    def __init__(self, **data):
        super().__init__(**data)

        self.structure_path = Path(self.structure_path) 
        # self.structure_path = self.structure_path / self.structure_database.lower()

        self.data_dir = Path(self.data_dir) if self.data_dir is not None else None
        
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
    """
    TODO: 
    - separate out a 'phosphositedataset' similar to used in CASM that can be used to compose a MotifAnalysis object
    - i.e. decouple a lot of the dataset stuff from MotifAnalysis and put in a separate class.
    """

    def __init__(
        self,
        config: MotifAnalysisConfig = None,
    ) -> None:

        self.dataset_path = config.dataset_path

        self.alphafold_structure_dir = config.alphafold_structure_dir if config.alphafold_structure_dir is not None else config.structure_path

        self.use_dataset = config.use_dataset

        self.af_model_params = config.af_model_params

        self.af_ignore_fragments    = config.af_ignore_fragments
        self.af_model_version       = config.af_model_version
        self.af_file_extension      = config.af_file_extension
 
        self.dbptm_path = config.dbptm_path

        self.include_isoforms = config.include_isoforms

        # Filtering
        self.species_filter = config.species_filter
        self._max_sites = config.max_sites

        self.radius = config.radius

        self.sites = []

        if self.use_dataset == "dbptm":
            self._load_dbptm()
        else:
            raise ValueError(f"Invalid dataset: {self.use_dataset}")

        # Load in structures 
        self._load_structures()

    
    def _load_dbptm(self) -> pd.DataFrame:
        """Load the DBPTM database."""
        
        if self.dataset_path is None: # use default dbPTM if `dataset_path` unspecified
            self.dataset_path = self.dbptm_path
        
        print(f"Loading from {self.dataset_path}")
        df = pd.read_csv(
            self.dataset_path, 
            sep="\t",
            names=[
                "entry_name",       # dbPTM columns
                "acc_id", 
                "position", 
                "mod_type",
                "pmids", 
                "seq_window",  
            ],
        )

        #print(df.entry_name)
        #exit(0)
        # Add species column, which is the 2nd part of `entry_name` 
        df["entry_name"] = df["entry_name"].astype(str)
        df["species"] = df["entry_name"].apply(lambda x: x.split("_")[-1].lower())

        # Filter by species
        if self.species_filter is not None:
            df = df[df["species"].isin(self.species_filter)]

        df = df.head(self._max_sites) # Restrict number of sites in dataset

        if not self.include_isoforms:
            # Remove rows where the `acc_id` is a seqeuence isoform.
            df = df[~df["acc_id"].str.contains("-")]

        
        # Create residue column from seq_window.  
        # This is the character in the `seq_window` which is calculated
        # by dividing the length of the `seq_window` by 2 and rounding down.
        df["seq_window"] = df["seq_window"].astype(str)
        df["residue"] = df["seq_window"].apply(lambda x: x[int(len(x)/2)]) 
        

        """
        TODO: 
        - for each entry, get sequence from uniprot `acc_id`
        - confirm that middle of `seq_window` residue matches the `acc_id` sequence at the `position` 
        
        """

    
        # Create list of sites for each PTM site recorded in the dataset. 
        for i, row in df.iterrows():
            ptm_site = PTMSite(
                entry_name=row["entry_name"],
                acc_id=row["acc_id"], 
                residue=row["residue"], 
                position=row["position"], 
                ptm_type=row["mod_type"], 
            )
            self.sites.append(ptm_site)

        self.dataset = df
        return df

    def _load_structures(self) -> None:

        # If use alphafold, load structures from alphafold directory

        # TODO: If use alternative structure database, load structures from that directory 
        
        failed_sites = []
        motifs = {}
        for site in tqdm(self.sites):
            #site.load_structure(self.alphafold_structure_dir, self.af_model_params)
            try:
                motif: StructuralMotif = self._load_alphafold(site=site)
                motifs[f"{site.entry_name}-{site.node_id}"] = motif
            except FileNotFoundError:
                failed_sites.append(site)
                print(f"Alphafold structure not found for {site.acc_id}")
                continue

        self.motifs = motifs
        self.failed_sites = failed_sites

    def _load_alphafold(
        self,
        site: PTMSite,

    ) -> None:

        """Load an alphafold structure."""

        

        #self.sites: List[PTMSite]          # set of sites (i.e. ID and residue position)
        #self.motifs: List[StructuralMotif]  # subgraph centred around each site # TODO make 
                                            # this efficient without multiple subgraphs; but 
                                            # rather just a refernce to the original graph 

        # TODO: maybe change from List to Dict with key as site ID for easier referencing? 

        
        pdb_path = self.alphafold_structure_dir / self._get_af_filename(site.acc_id)
        print(pdb_path)

        # Check if path is an existing file
        if not pdb_path.is_file():
 
            raise FileNotFoundError(f"Alphafold structure not found for {site.acc_id} at {pdb_path}")
        
        motif = StructuralMotif(
            site=site,
            structure_path=pdb_path, 
            radius=self.radius,
        )
        return motif
            

        



    def _get_af_filename(
        self,
        acc_id: str, 
    ) -> str: 

        if self.af_ignore_fragments:
            fragment = 1
            return f"AF-{acc_id}-F{fragment}-model_v{self.af_model_version}.{self.af_file_extension}"
        else:
            raise NotImplementedError(f"Multiple fragments for AF structure not implemented.")
    

    """
    Analysis methods
    """
    def difference_transform(
        self,
    ) -> None:
        """Perform a difference transform on the set of motifs."""

        pass


