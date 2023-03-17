"""Base config object for use with Motif analysis."""

from pydantic import BaseModel


class MotifGraphConfig(BaseModel):
    """Config for motif graph analysis."""
    
    granularity: str = "residue"
