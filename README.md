# hopa - Homogenized orientation using pairwise alignments

## Rationale
De novo sequence assembly and structural variant analyses can generate lists of sequences that are similar, but in arbitrary orientations. HOPA simply takes an input FASTA file and aligns each sequence (forward and reverse complement) against the first sequence in the list to orient all respectively for downstream analyses.
