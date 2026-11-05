#!/usr/bin/env bash
set -euo pipefail

# Split per-gene alignments into backbone and query fasta files (unaligned)
# $1 alignment dir (per-gene fasta files)
# $2 backbone tree
# $3 output dir

ALNDIR=$1
BBONE=$2
OUTDIR=$3

BBDIR=$OUTDIR/placement/backbone
QDIR=$OUTDIR/placement/query
mkdir -p $BBDIR $QDIR

# get backbone taxa from tree
nw_labels -I $BBONE | sort -u > $OUTDIR/placement/backbone_taxa.txt

for aln in "$ALNDIR"/*; do
    gene=$(basename "$aln")

    # split into backbone and query per gene
    seqkit grep -f $OUTDIR/placement/backbone_taxa.txt -w 0 --quiet "$aln" -o $BBDIR/$gene
    seqkit grep -v -f $OUTDIR/placement/backbone_taxa.txt -w 0 --quiet "$aln" -o $QDIR/$gene
done

# backbone tree: strip internal labels, resolve polytomies
nw_topology -bI $BBONE > $OUTDIR/placement/backbone_stripped.tree
python -c "import treeswift as ts; t=ts.read_tree_newick('$OUTDIR/placement/backbone_stripped.tree'); \
          [c.resolve_polytomies() for c in t.root.children]; print(t)" > $OUTDIR/placement/backbone.tree
rm $OUTDIR/placement/backbone_stripped.tree
