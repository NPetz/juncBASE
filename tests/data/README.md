# Real-data test fixtures

The GTF records are exons from FlyBase *Drosophila melanogaster* release
6.69 (`dmel-all-r6.69.gtf`), for roX2 (FBgn0019660) and dre4 (FBgn0002183).
The two BAM files are extracts from the public Pasilla RNA-seq runs SRR031716
and SRR031708, aligned by the associated Pasilla project against the same
reference. Each includes alignments overlapping `3L:1871574-1876312` and
`X:11579065-11580100` (1-based, inclusive). SRR031716 is paired-end and
SRR031708 is single-end. The names and mate flags are preserved; the extracts
are not intended as complete read-pair datasets.

The cross-version test supplies empty intron-boundary inputs, so it is a
software regression check, not evidence of differential intron retention or
correct fragment-level counting. The GTF contig names and attribute order are
transformed in the test so both versions receive identical, upstream-compatible
inputs; these changes do not alter the original coordinates.
