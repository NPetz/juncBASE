# Python 3 port: work in progress

This branch is an experimental port, not a validated replacement for JuncBASE.
The upstream master branch remains the reference for existing behavior.

The main scripts and their coordinate-tree dependencies have received an
initial Python 3 conversion, using Python 3.11's `lib2to3` and targeted
compatibility changes. This is not yet a complete or validated port.
`clusterASExons2.py` could not be automatically converted because its upstream
source contains an incomplete expression. Optional dependencies and some
legacy wrapper scripts also remain unresolved.

One smoke test used Python 3.13 and pysam 0.24.1 on 102 alignments extracted
from the real Pasilla run SRR031716 around `X:11579065-11580100`. The command
completed and emitted genome-read and junction files; two reported junctions
had counts of 2 and 3, consistent with the unclipped CIGAR records inspected.
For this extract, the original Python 2 and ported preprocessors produced
identical junction tables and identical decompressed genome-read data.
This does not establish full-pipeline equivalence.

The Python 3 event caller reached completion in an exploratory software smoke
test using the same sample twice and empty intron-boundary inputs. The annotation
builder interpreted gene symbols as transcript IDs for the tested FlyBase GTF,
and the BAM and annotation chromosome names did not agree; therefore this run
does not constitute a valid biological analysis or a numerical regression test.

Important unresolved issues include paired-end fragment counting, chromosome
prefixes, GTF attribute parsing, clipped/indel alignments, and compatibility
of remaining dependencies. Do not use this branch for scientific conclusions
until these issues and full regression comparisons have been addressed.
