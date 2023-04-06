"""Classes for motifs."""

import networkx as nx
import click as ck 

from kimono.motif.definitions import DEFAULT_PROTEIN_GRAPH_CONFIG
from kimono.ptm import PTMSite

from pathlib import Path

from graphein.protein import ProteinGraphConfig 
from graphein.protein.graphs import construct_graph 
from graphein.protein.config import DSSPConfig

from graphein.protein.subgraphs import extract_subgraph_from_point
from graphein.protein.edges.distance import *	
from graphein.protein.features.nodes import rsa as rsa_function

class LinearMotif():
    """Represents a sequence of residues that form a motif.  The actual
    position of residues in 3D structure are not meaningful.
    """
    pass 

class PsuedoLinearMotif():
    """Represents a 3D structural motif that is analagous to 
    a linear motif in terms of structural pattern, but is not 
    composed of consecutive sequence-adjacent residues. 
    """
    pass 

class StructuralMotif():
    """Represents a structural motif on a protein.
    
    It is a graph of residues (or atomistic) that are within a threshold distance from the 
    centre node, which is a post-translational modification (PTM) site.

    The 'bubble' around the centre node is the motif.  The radius can be updated. 
    """

    """
    TODO: 
    - we need a directory to store these motifs (serialisable) so that we don't have
    - multiple copies of the same graph (protein ID) in memory. 
    - each instance of this class that references the same protein ID should point to the same graph; 
    and simply return a different subgraph. 
    - obviously the first time it needs to create the graph. 
    """

    def __init__(
        self, 
        g: nx.Graph = None, 
        site: PTMSite = None,
        structure_path: Path = None,
        filename: str = None,
        graph_config: ProteinGraphConfig = DEFAULT_PROTEIN_GRAPH_CONFIG,
        radius: float = 12.0, # Distance threshold from the center node in Ångströms 
        rsa: float = 0.0, # Relative solvent accessibility threshold 
        granularity: str = "residue", # TODO: implement atomistic granularity
    ) -> None:
        
        if site is None:
            raise ValueError("Must provide a PTM site.")

        self.centre_node = site.node_id 
        self.granularity = granularity

        self._radius    = radius
        self._rsa       = rsa

        if g is None:
            # TODO: create graph from uniprot ID and residue position 
            if structure_path is None:
                raise ValueError("Must provide either a graph or a structure path.")
            
            
            # Check if structure is directory
            if structure_path.is_dir():
                if filename is None:
                    raise ValueError("Must provide a filename for structure path.")
                structure_path = structure_path / filename 

            # Check if structure_path exists 
            if not structure_path.exists():
                raise ValueError(f"Path does not exist: {structure_path}") 

            g = construct_graph(pdb_path=structure_path, config=graph_config)

        self.g = g 

        self._motif = self._get_motif_subgraph() # GET SUBGRAPH 

        

    @property
    def radius(self) -> float:
        return self._radius

    @radius.setter
    def radius(self, radius: float) -> None:
        self._radius = radius

        # Reload subgraph
        self._motif = self._get_motif_subgraph()

        """
        TODO:
        - reload subgraph from original graph 
        """

    @property
    def motif(self) -> nx.Graph:
        return self._motif

    @property
    def nodes(self):
        return self._motif.nodes(data=False)

    def _get_motif_subgraph(self) -> nx.Graph:
        """Get the subgraph of the motif."""
        try:
            g_site: Dict = self.g.nodes(data=True)[self.centre_node]
        except: 
            raise ValueError(f"Centre node '{self.centre_node}' not found in graph.")
        
        g_site_rsd: str = str(g_site['residue_name'])
        g_site_pos: int = int(g_site['residue_number'])

        # Subgraph (radius)
        s_g = self._get_protein_subgraph_radius(g=self.g, site=self.centre_node, r=self.radius)

        # Subgraph (rsa)
        # TODO

        return s_g

    @property
    def nodes(self):
        return self._motif.nodes(data=False)

    def _get_protein_subgraph_radius(self, g: nx.Graph, site: str, r: float) -> nx.Graph:
        """Get the subgraph of the motif."""
        # get centre point   
        try:
            x_y_z = node_coords(g, site)
        except ValueError:
            raise ValueError("Specified site isn't in correct format.")       
        
        # Get subgraph
        s_g = extract_subgraph_from_point(g, centre_point=x_y_z, radius=r)
        return s_g

    def average_difference_transform(self):
        """Average difference transform the motif."""
        
        diff = self.difference_transform()
        return sum(diff) / len(diff)
        
    def sum_difference_transform(self):
        """Sum difference transform the motif."""
        
        diff = self.difference_transform()
        return sum(diff)

    def max_difference_transform(self):
        """Max difference transform the motif."""
        
        diff = self.difference_transform()
        return max(diff)

    def difference_transform(
        self,
        zeroed: bool = True, # Consecutive residues will result in 0 difference 
    ):
        """Difference transform the motif."""
        
        # Iterate through all nodes in motif 
        l = [
            int(n.split(":")[-1])
            for n in self.nodes
        ]
        l.sort()

        # Calculate difference transform
        return [l[i] - l[i-1] - int(zeroed) for i in range(1, len(l))]

    def pos_difference_transform(
        self,
    ):
        """Difference transform with only non zero values."""
        diff = self.difference_transform()
        return [d for d in diff if d != 0]
            
    def __repr__(self) -> str:
        return f"StructuralMotif({self.g.name} @ {self.centre_node}, granularity={self.granularity}, radius={self.radius})"


    
        
