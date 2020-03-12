# SPeQ - Strain PErcent Quantifier

## Rationale
Metagenomic sequencing studies can contain DNA from multiple species, which are often relatively easy to distinguish. However, if a sample contains multiple variants of the same virus/organism, a large fraction of the reads may be synonymous in both variants. SPeQ only uses the polymorphisms that distinguish individual variants in order to sensitively quantify percentages of strains/variants.


## Installation
```
git clone --recursive https://github.com/ibharvey/speq
cd speq
mkdir -p build
cd build
cmake ..
#Builds a binary located at ..../speq/build/src/speq
make
#Installs the binary to be used directly from the commandline
#   (i.e. 'speq [options]' instead of ./speq [options]
sudo make install

```

SPeQ requires Seqan3 and its respective dependencies. FYI this includes a fairly new compiler. I have only successfully worked with Seqan3 on Ubuntu 18

## Algorithm

TotalNumberOfKmers * Prob(NoSequencingErrorsPerKmer) * **VariantPercentage** * PercentageOfUniqueKmers = NumberOfUniqueKmers

- TotalReads is known for a given FAST(AQ) file.
- PercentUniqueKmers can be found by iterating a sliding kmer window over each reference sequence, back-checking against all reference sequences.
- NumberUniqueKmers can be found by iterating a sliding kmer window through all reads, only counting kmers that exclusively match one variant.

