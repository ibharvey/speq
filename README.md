# speq - Strain PErcent Quantifier

## Rationale
Metagenomic sequencing studies can contain DNA from multiple species, which are often relatively easy to distinguish. However, if a sample contains multiple variants of the same virus/organism, a large fraction of the reads may be synonymous in both variants. SPeQ only uses the polymorphisms that distinguish individual variants in order to sensitively quantify percentages of strains/variants.


## Installation
```
git clone --recursive https://github.com/ibharvey/speq
cd speq
mkdir -p build
cd build
cmake ..
make
#Builds a binary located at ..../speq/build/src/speq
```

SPeQ requires Seqan3 and its respective dependencies. FYI this includes a fairly new compiler.
