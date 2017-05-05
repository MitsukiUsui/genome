# What's this?
This is a sub project of altorf, analyzing altorf from genome perspective. It creates simulated genomes for each species and compares altorf properties with real genome. Programs are organized into further sub-directories depending on the stage of analysis. Basically written in C++ or Python3, and also Shell Scripts are used to create pipelines.

## Workflow
Follow this flow from top to bottom. For further information, please refer to README on each sub directories.

### speciespick
* Chooses which species in RefSeq to analyze.

### preprocess
* Splits downloaded .fna files and fills ambiguous bases randomly if needed.

### genomeshuffle
* Creates simulated genomes according to .fna & .gff files.

### patternanalyze
* Analyze pattern of altorf comparing real genome and simulated genome.


## Prerequisites
I chose SeqAn as C++ library in order to process DNA sequence efficiently and to make biological fileI/0 easier.
* seqan: <https://www.seqan.de/>
    * An open source C++ library of efficient algorithms and data structures for the analysis of sequences.
    * For installation, please refered to <http://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html>  
        0. Mac: `brew install homebrew/science/seqan`
        0. servers: clone github repository <https://github.com/seqan/seqan>
* CMake: <https://cmake.org/>
    * seqan recommends compilation with CMake <http://seqan.readthedocs.io/en/master/Infrastructure/Use/FindSeqAnCMake.html>

* Anaconda (ver 3.X)
* BioPython

## Compilation of .cpp
* Commands below creates binary file in build directory.
```
cd $PATH_TO_DESIGNATED_DIR
mkdir build
cd build
cmake ../src
make
```
* In case you obtained SeqAn from a git clone, you need to specify the install location and include path.
    * `cmake ../src -DCMAKE_PREFIX_PATH="$HOME/software/seqan/util/cmake" -DSEQAN_INCLUDE_PATH="$HOME/software/seqan/include"`
        * Please change paths according to your environment.
