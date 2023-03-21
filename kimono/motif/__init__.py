"""Classes for motifs."""

import networkx as nx
import click as ck 


class StructuralMotif():

    def __init__(
        self, 
        g: nx.Graph = None,

    ) -> None:

        self.g = g if g is not None else nx.Graph()
    
    
