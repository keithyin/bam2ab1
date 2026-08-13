#!/bin/bash
# Regenerate test-data/empty-window/ (seq.fasta, query.fasta, empty-window.sort.bam).
#
# The reference is 1200bp (contig "seq1"); the 4 query reads (300bp, exact
# substrings) align onto [0,600) only. With `--chunkSize 300`, windows
# [600,900) and [900,1200) have zero read coverage, which exercises the
# empty-window branch in plp_from_records (the F2 fix, see lines 55-68).
# The same data is asserted by test_plp_from_records_empty_window_from_bam.
#
# Usage:
#   MINIMAP2=/path/to/minimap2 ./scripts/make_empty_window_test.sh
# (defaults to `minimap2` on PATH; if you lack root, grab the binary with:
#   apt-get download minimap2 && dpkg-deb -x minimap2_*.deb . \
#     && MINIMAP2=$PWD/root/usr/bin/minimap2 ./scripts/make_empty_window_test.sh)

set -euo pipefail

DIR="$(cd "$(dirname "$0")/.." && pwd)/test-data/empty-window"
MINIMAP2="${MINIMAP2:-minimap2}"
mkdir -p "$DIR"
cd "$DIR"

python3 - <<'PY'
import random

random.seed(42)
seq = ''.join(random.choice('ACGT') for _ in range(1200))

def wrap(s, width=80):
    return '\n'.join(s[i:i + width] for i in range(0, len(s), width))

with open('seq.fasta', 'w') as f:
    f.write(f'>seq1\n{wrap(seq)}\n')

with open('query.fasta', 'w') as f:
    for idx, off in enumerate([0, 100, 200, 300], start=1):
        f.write(f'>read{idx}\n{wrap(seq[off:off + 300])}\n')
PY

"$MINIMAP2" -ax sr --secondary=no seq.fasta query.fasta 2>/dev/null \
    | samtools sort -@ 4 -o empty-window.sort.bam -
samtools index empty-window.sort.bam

echo "Wrote: $DIR/{seq.fasta, query.fasta, empty-window.sort.bam, empty-window.sort.bam.bai}"
