# Kimono

> **ki**nase **mo**tif **no**tation

## First pass 

- Scan all 3D structures of phosphosites (pLDDT > threshold)
- Compare PDB structures with AF where structure exist 
- first, generate a sequence for each *N* Å bubble, where the sequence numbers start from 0 (each residue).  The sequential nature of these can be visualised (i.e. is it 0, 1, 2, 3, 4 or 0, 1, 10, 11 etc. ) 
- either order by `RES` seq position, or by increasing distance from `MOD_RSD` 
- see if this correlates with pLDDT of (psite | average in bubble) which i suspect it will

- similar to `structuremap`, develop a score for each site along with its `in.IDR` etc. scores 
- this score represents the "non-linearity" of the particular motif 

The score should reflect the following:

1. if all the residues within the bubble are sequence-continuous (i.e. 12345...), the score should be 0. 
1. if there are adjacent residues in space (euclidean distance) that are far apart in sequence, the score should be high.  
1. if the 



Use difference transform 
