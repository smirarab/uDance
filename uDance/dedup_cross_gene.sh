#!/usr/bin/env bash
set -euo pipefail

# Cross-gene deduplication: identify duplicates only, no fasta rewriting.
# Outputs dedup_map.txt and per-gene dupmap/remove files.
# Downstream split_bb_query_per_gene.sh uses these to exclude duplicates.
#
# $1 input alignment dir (per-gene fasta files)
# $2 output metadata dir (dedup_map.txt, per-gene .dupmap.txt/.remove.txt)
# $3 backbone tree (to extract backbone taxa list)

INDIR=$1
OUTDIR=$2
BBONE=$3

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

mkdir -p "$OUTDIR"

# extract backbone taxa from tree
nw_labels -I "$BBONE" | sort -u > "$TMPDIR/backbone_taxa.txt"

# Phase 1: identify per-gene duplicates (read-only)
echo "Phase 1: identifying per-gene duplicates..."
for fa in "$INDIR"/*; do
    [ -f "$fa" ] || continue
    gene=$(basename "$fa")
    seqkit rmdup -s -D "$TMPDIR/$gene.dup.txt" "$fa" > /dev/null 2>&1 || true
done

# Phase 2: compile cross-gene duplicates and write metadata
echo "Phase 2: compiling cross-gene duplicates..."
python -I "$SCRIPT_DIR/compile_cross_gene_dups.py" "$TMPDIR" "$OUTDIR" "$INDIR" "$TMPDIR/backbone_taxa.txt"

echo "Cross-gene deduplication identification complete. Metadata in $OUTDIR"
