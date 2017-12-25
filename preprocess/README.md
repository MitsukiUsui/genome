# preprocess

## Summary
Splits downloaded .fna file into two files (chromosome and plasmid).  
Fill ambiguous bases if needed.  
For further information on usage, please refer to split_fasta.html or `./split_fasta -h`  

## Input
* Downloaded .fna file from NCBI RefSeq

## Scripts
1. `./scripts/rename.py`
    * gunzip downloaded `.fna`, `.gff`, and `.cds`
    * add alias using taxid as a unique identifier

1. `./scripts/caller.sh`
    * divide downloaded `.fna` into chromosomal or plasmid sequences.
    * fill Ns randomly
