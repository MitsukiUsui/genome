# genomeshuffle

## Summary
Creates simulated genome according to `.fna` & `.gff`  
For further information on usage, please refer to genome_shuffle.html or `./genome_shuffle -h`

## Input
* Preprocessed `.fna`
* Downloaded `.gff` from RefSeq

## Compilation of .cpp
* Commands below creates binary file in build directory.
```
cd $PATH_TO_DESIGNATED_DIR
cd build
cmake ../ 
make
```
* If you fail setting CMAKE_CXX_COMPILER properly, recreate build directory. It can only be set the first time cmake is run in a given build directory.
* In case you obtained SeqAn from a git clone, you need to specify the install location and include path.
    * `cmake ../ -DCMAKE_PREFIX_PATH="$HOME/software/seqan/util/cmake" -DSEQAN_INCLUDE_PATH="$HOME/software/seqan/include"`
        * Please change paths according to your environment.


## Implementation Details

### Class
1. GeneticCode
    * A class to count occurance of each codon in order to achieve synonymous shuffling according to codon usage
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
    * get a vector of ShuffleRegion according to　myCDS_vecvec
1. shuffle_genome()
    * shuffle sequences according to shuffleRegions return by get_shuffle_region()
    * 3 shuffing method (base, codon, synonymous) is defined so far
