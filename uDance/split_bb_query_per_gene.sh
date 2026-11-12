#!/usr/bin/env bash
set -euo pipefail

# Split per-gene fastas into backbone and query, filtering out cross-gene duplicates.
# $1 alignment dir (original per-gene fasta files)
# $2 backbone tree
# $3 output dir
# $4 dedup metadata dir (optional; contains dedup_map.txt, per-gene .remove.txt)

ALNDIR=$1
BBONE=$2
OUTDIR=$3
DEDUPDIR=${4:-}

BBDIR=$OUTDIR/placement/backbone
QDIR=$OUTDIR/placement/query
mkdir -p $BBDIR $QDIR

# get backbone taxa from tree
nw_labels -I $BBONE | sort -u > $OUTDIR/placement/backbone_taxa.txt

# if dedup metadata exists, build global remove list and prune backbone taxa list
if [ -n "$DEDUPDIR" ] && [ -f "$DEDUPDIR/dedup_map.txt" ]; then
    cut -f2- "$DEDUPDIR/dedup_map.txt" | tr '\t' '\n' | sort -u > "$OUTDIR/placement/dedup_remove.txt"
    grep -vxFf "$OUTDIR/placement/dedup_remove.txt" "$OUTDIR/placement/backbone_taxa.txt" \
        > "$OUTDIR/placement/backbone_taxa_deduped.txt" || true
    mv "$OUTDIR/placement/backbone_taxa_deduped.txt" "$OUTDIR/placement/backbone_taxa.txt"
    # copy dedup_map.txt for downstream use (expand_dedupe_newick)
    cp "$DEDUPDIR/dedup_map.txt" "$OUTDIR/placement/dedup_map.txt"
fi

for fa in "$ALNDIR"/*; do
    [ -f "$fa" ] || continue
    gene=$(basename "$fa")

    # per-gene remove list from dedup
    remove_file="$DEDUPDIR/$gene.remove.txt"

    if [ -n "$DEDUPDIR" ] && [ -f "$remove_file" ] && [ -s "$remove_file" ]; then
        # backbone: grep backbone taxa, then exclude duplicates
        seqkit grep -f $OUTDIR/placement/backbone_taxa.txt -w 0 --quiet "$fa" \
            | seqkit grep -v -f "$remove_file" -w 0 --quiet -o $BBDIR/$gene
        # query: grep non-backbone, then exclude duplicates
        seqkit grep -v -f $OUTDIR/placement/backbone_taxa.txt -w 0 --quiet "$fa" \
            | seqkit grep -v -f "$remove_file" -w 0 --quiet -o $QDIR/$gene
    else
        seqkit grep -f $OUTDIR/placement/backbone_taxa.txt -w 0 --quiet "$fa" -o $BBDIR/$gene
        seqkit grep -v -f $OUTDIR/placement/backbone_taxa.txt -w 0 --quiet "$fa" -o $QDIR/$gene
    fi
done

# backbone tree: strip internal labels, prune deduped taxa, resolve polytomies
nw_topology -bI $BBONE > $OUTDIR/placement/backbone_stripped.tree

if [ -f "$OUTDIR/placement/dedup_remove.txt" ] && [ -s "$OUTDIR/placement/dedup_remove.txt" ]; then
    nw_prune $OUTDIR/placement/backbone_stripped.tree $(cat "$OUTDIR/placement/dedup_remove.txt") \
        > $OUTDIR/placement/backbone_pruned.tree
    mv $OUTDIR/placement/backbone_pruned.tree $OUTDIR/placement/backbone_stripped.tree
fi

python -c "import treeswift as ts; t=ts.read_tree_newick('$OUTDIR/placement/backbone_stripped.tree'); \
          [c.resolve_polytomies() for c in t.root.children]; print(t)" > $OUTDIR/placement/backbone.tree
rm $OUTDIR/placement/backbone_stripped.tree
