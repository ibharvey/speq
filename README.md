# SPeQ - Strain PErcent Quantifier

## Rationale
Metagenomic sequencing studies can contain DNA from multiple species, which are often relatively easy to distinguish. However, if a sample contains multiple variants of the same virus/organism, a large fraction of the reads may be synonymous in both variants. SPeQ only uses the polymorphisms that distinguish individual variants in order to sensitively quantify percentages of strains/variants. Because of this, it can be extremely slow compared to other methods, but maintains high sensitivity.

## Usage
### SPeQ-Index
Organize a directory so that each of its subdirectories is a genomic variant, and each FASTA/Q file in that subdirectory represents one isolate's sequence. 


### SPeQ-Scan
Compare reads to defined genotypes to determine the fractional abundance of each in a sample.

## Installation
### With Docker (from source)
If you are using a Linux/Mac distribution and have a Docker daemon running, you can install with:
```
git clone --recursive https://github.com/ibharvey/speq
bash speq/dockerfiles/alpine/setup.sh speq/
```
The benefit being that docker will pull all of the necessary dependencies into a separate container from your OS, giving you less of a headache.

To run speq through this Docker image, use:
```
docker run speq [index|scan|all] [options]
```

### With CMake (from source)
```
git clone --recursive https://github.com/ibharvey/speq
cd speq
mkdir -p build
cd build
cmake ..
make
sudo make install

```
SPeQ requires Seqan3 and its respective dependencies. FYI this includes a fairly new compiler. 

To run speq following CMake installation, use:
```
speq [index|scan|all] [options]
```

## Algorithm

TotalNumberOfKmers * Prob(NoSequencingErrorsPerKmer) = NumberOfPerfectKmers

NumberOfPerfectKmers * **VariantPercentage** = NumberOfVariantSpecificKmers

NumberOfVariantSpecificKmers * PercentageOfUniqueKmersInVariant = NumberOfVariantUniqueKmers

- TotalReads is known for a given FAST(AQ) file.
- PercentUniqueKmers can be found by iterating a sliding kmer window over each reference sequence, back-checking against all reference sequences.
- NumberUniqueKmers can be found by iterating a sliding kmer window through all reads, only counting kmers that exclusively match one variant.

