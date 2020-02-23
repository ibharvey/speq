# hopa - Homogenized orientation using pairwise alignments

## Rationale
De novo sequence assembly and structural variant analyses can generate lists of sequences that are similar, but in arbitrary orientations. HOPA simply takes an input FASTA file and aligns each sequence (forward and reverse complement) against the first sequence in the list to determine the 'best' orientation of each sequence to the first and outputs this oriented file. 



## Installation
```
git clone --recursive https://github.com/ibharvey/hopa
cd hopa
mkdir -p build
cd build
cmake ..
make
#Builds a binary located at ..../hopa/build/src/hopa
```

HOPA requires Seqan3 and its respective dependencies. FYI this includes a fairly new compiler.
