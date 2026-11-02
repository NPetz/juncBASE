# Python 3 port: work in progress

This branch is an experimental port, not a validated replacement for JuncBASE.
The upstream master branch remains the reference for existing behavior.

The main scripts and their coordinate-tree dependencies have received an
initial Python 3 conversion, using Python 3.11's `lib2to3` and targeted
compatibility changes. This is not yet a complete or validated port.
`clusterASExons2.py` contained two upstream syntax errors and expected an
older event-table format. With approval, it now reuses the current count-table
parser; its raw and normalized output matches the current table builder on the
real-data test below. The normalized mode now requires `--jcn_seq_len`.
Other optional dependencies, legacy wrapper scripts and unused vendored `opus7`
modules remain unresolved.

One smoke test used Python 3.13 and pysam 0.24.1 on 102 alignments extracted
from the real Pasilla run SRR031716 around `X:11579065-11580100`. The command
completed and emitted genome-read and junction files; two reported junctions
had counts of 2 and 3, consistent with the unclipped CIGAR records inspected.
For this extract, the original Python 2 and ported preprocessors produced
identical junction tables and identical decompressed genome-read data.

The opt-in `tests/test_python2_equivalence.py` compares real BAM extracts at
roX2 and dre4 (SRR031716 paired-end and SRR031708 single-end) and FlyBase
r6.69 exon records. It checks annotation rows, junction and genome-read
outputs, exon–intron boundary counts, event records and counts (cassette exons
and intron retention), Fisher-test output, and a non-empty aggregated count
table. It also checks the
updated cluster script against raw and normalized current-format count tables.
The aggregation
smoke test duplicates the same generated event table under two sample labels;
it is not a comparison of independent biological replicates.
The test normalizes the annotation's
chromosome prefix and attribute order for the original parser. Because the old
`pysqlite2` package could not be built in the reference environment, its
temporary Python 2 checkout uses Python 2's `sqlite3` in place of that import.
The boundary counts are generated from the same real BAM extracts. This
comparison is not a differential splicing analysis or a full validation of
retention counting, paired-end fragment handling, or all event types.
Ordering of event and junction records is ignored.

To run the cross-version test, set `JUNCBASE_PY2_ROOT` and
`JUNCBASE_PY2_PYTHON` to an independently prepared reference checkout and
Python 2 environment, respectively, and run `python -m unittest discover -s
tests -v` in a Python 3 environment with Biopython, pysam, R and rpy2.
Each environment must contain its own working R and rpy2 installation.

The Python 3 event caller reached completion in an exploratory software smoke
test using the same sample twice and empty intron-boundary inputs. At the time,
the annotation builder interpreted gene symbols as transcript IDs and the BAM
and annotation chromosome names did not agree. This run does not constitute a
valid biological analysis or a numerical regression test.

The annotation builder now retrieves `transcript_id` by name rather than by
attribute position. An integration test with real FlyBase r6.69 roX2 exon records
checks that two transcripts and their coordinates are stored separately. This
fix changes the legacy behavior for GTFs whose attributes are ordered differently.

Subprocess failures in `getASEventReadCounts.py` now stop the analysis instead
of silently producing incomplete outputs. This is an intentional change from
upstream error handling; successful runs retain the same counting logic.

Important unresolved issues include paired-end fragment counting, chromosome
prefixes, other GTF parsing assumptions, clipped/indel alignments, and compatibility
of remaining dependencies. Do not use this branch for scientific conclusions
until these issues and full regression comparisons have been addressed.
