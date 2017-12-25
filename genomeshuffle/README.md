# genomeshuffle

## Summary
Creates simulated genome according to `.fna` & `.gff`  
For further information on usage, please refer to genome_shuffle.html or `./genome_shuffle -h`

## Input
* Preprocessed `.fna`
* Downloaded `.gff` from RefSeq

## Implementation Details

### Class
1. GeneticCode
    * A class to achieve synonymous shuffling according to codon usage
1. MyCDS
    * A wrapper class for seqan GffRecord to store information about the cds (label and phase)
1. ShuffleRegion
    * A container which stores nformation needed for shuffling.
1. ShuffleRegionFactory
    * A factory class which create vector of ShuffleRegion according to 

### method
1. get_myCDS_vecvec()
   * get a vector of vector of MyCDS, each vector contains MyCDSs which belongs to each sequence
1. get_shuffle_region()
    * get a vector of ShuffleRegion according myCDS_vecvec
1. shuffle_genome()
    * shuffle sequences according to shuffleRegions return by get_shuffle_region()
    * 3 shuffing method is defined so far
        * shuffle_base()
        * shuffle_codon()
        * and shuffle_synonymous()
