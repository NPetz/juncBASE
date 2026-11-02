# Python 3 port: work in progress

This branch is an experimental port, not a validated replacement for JuncBASE.
The upstream master branch remains the reference for existing behavior.

Only `preProcess_getASEventReadCounts.py` has been converted so far, using
Python 3.11's standard `lib2to3` tool. Its gzip output is opened in text mode
for Python 3. The rest of the event-calling pipeline still needs porting.

One smoke test used Python 3.13 and pysam 0.24.1 on 102 alignments extracted
from the real Pasilla run SRR031716 around `X:11579065-11580100`. The command
completed and emitted genome-read and junction files; two reported junctions
had counts of 2 and 3, consistent with the unclipped CIGAR records inspected.
This does not establish equivalence to the original Python 2 implementation.

Important unresolved issues include paired-end fragment counting, chromosome
prefixes, GTF attribute parsing, clipped/indel alignments, and compatibility
of remaining dependencies. Do not use this branch for scientific conclusions
until these issues and full regression comparisons have been addressed.
