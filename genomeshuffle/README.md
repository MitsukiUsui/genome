# genomeshuffle

## Summary
create simulated genome according to .fna & .gff
For further information on usage, please refer to genome_shuffle.html or `./genome_shuffle -h`

## Input
* preprocessed .fna
* downloaded .gff from RefSeq

## Implementation Details

### GeneticCode
#### Summary
* Class to achieve synonymous shuffling according to codon usage.
#### DataStructure
* codons and amino acids are each managed by codonIds (0 ~ 63) and aaIds (0~22)
    * aaId0 corresponds to stop

0. codonToAa[64]
    * given codonId, return corresponding amino acid
0. aaIndex[22] & aaToCodon[64]
    * given aaId, the location of corresponding codons are given by [aaIndex[aaId], aaIndex[aaId + 1])
        * therefore,  this is always the case that `aaIndex[0]==0 & aaIndex[22]==64`
    * you can specify corresponding codonIds by accessing aaToCodon with given range of index
* codonCount[64] & codonFreq[64]

### CDSs
#### Summary
