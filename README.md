# Genome
This is a sub project of altorf, focusing on genome.

## Summary
Creates Simulated genomes and compare altorf properties with real genome.

## Role of each Directory
0. speciespick: choose which species in RefSeq to analyze.
0. preprocess:  split downloaded .fna file, and fill ambiguous bases if needed
0. genomeshuffle: create simulated genome according to .gff file
0. patternanalyze: analyze pattern of altorf comparing real genome and simulated genome


## Prerequisites
* seqan: <https://www.seqan.de/>
    * an open source C++ library of efficient algorithms and data structures for the analysis of sequences
    * for installation, please refered to <http://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html>  
    in short:
        0. Mac...`brew install homebrew/science/seqan`
        0. Beyond...clone github repository
* CMake: <https://cmake.org/>
    * seqan recommends compilation with CMake <http://seqan.readthedocs.io/en/master/Infrastructure/Use/FindSeqAnCMake.html>

## Compilation off .cpp
0. create build directory and run `cmake ../src` in build
    * In case you obtained SeqAn from a git clone, you need to specify the install location and include path.
        * `cmake ../src -DCMAKE_PREFIX_PATH="$HOME/software/seqan/util/cmake" -DSEQAN_INCLUDE_PATH="$HOME/software/seqan/include"`
            * you need to change path according to your
0. `make` will compile designated .cpp and create binary in build
