# SPeQ - Strain PErcent Quantifier

## Rationale
Metagenomic sequencing studies can contain DNA from multiple species, which are often relatively easy to distinguish. However, if a sample contains multiple variants of the same virus/organism, a large fraction of the reads may be synonymous in both variants. SPeQ only uses the polymorphisms that distinguish individual variants in order to sensitively quantify percentages of strains/variants.


## Installation
### With Docker
If you are using a Linux/Mac distribution and have a Docker daemon running, you can install with:
```
git clone --recursive https://github.com/ibharvey/speq
cd speq
docker build . -t speq
```
The benefit being that docker will pull all of the necessary dependencies into the separate container, giving you less of a headache.

To run speq through the Docker image, use:
```
docker run speq [options]
```

### With CMake
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
SPeQ requires Seqan3 and its respective dependencies. FYI this includes a fairly new compiler. I have only successfully worked with Seqan3 on Ubuntu 18.04.4.

To run speq following CMake installation, use either:
```
speq [options]
#or
./speq [options]
```
Depending on how you installed the program.

## Algorithm

TotalNumberOfKmers * Prob(NoSequencingErrorsPerKmer) = NumberOfPerfectKmers

NumberOfPerfectKmers * **VariantPercentage** = NumberOfVariantSpecificKmers

NumberOfVariantSpecificKmers * PercentageOfUniqueKmersInVariant = NumberOfVariantUniqueKmers

- TotalReads is known for a given FAST(AQ) file.
- PercentUniqueKmers can be found by iterating a sliding kmer window over each reference sequence, back-checking against all reference sequences.
- NumberUniqueKmers can be found by iterating a sliding kmer window through all reads, only counting kmers that exclusively match one variant.

