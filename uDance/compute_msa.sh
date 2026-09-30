#!/usr/bin/env bash
set -euo pipefail

INP=$1      # raw unaligned per-gene fasta
OUT=$2      # destination: alignments/<gene> (what trimtaper expects)
TOOL=$3     # witch | twilight
CHAR=$4     # nuc | prot (uDance's existing chartype values)
THR=$5      # cpu threads
WORKDIR=$6  # scratch dir, unique per gene

mkdir -p "$WORKDIR"

case "$TOOL" in
  witch)
    molecule="dna"
    [ "$CHAR" == "prot" ] && molecule="amino"
    witch-msa -i "$INP" -d "$WORKDIR" -o aligned.fasta -t "$THR" --molecule "$molecule"
    mv "$WORKDIR/aligned.fasta" "$OUT"
    ;;
  twilight)
    type="n"
    [ "$CHAR" == "prot" ] && type="p"
    # Step 1: quick guide tree via MAFFT's --treeout side-output.
    mafft --retree 0 --treeout --reorder --thread "$THR" "$INP" > "$WORKDIR/mafft.out"
    python3 uDance/mafft2nwk.py "$INP.tree" "$INP" "$WORKDIR/guide.nwk"
    rm -f "$INP.tree" "$WORKDIR/mafft.out"
    # Step 2: TWILIGHT does the actual progressive alignment along that guide tree.
    twilight -i "$INP" -t "$WORKDIR/guide.nwk" -o "$OUT" -C "$THR" --type "$type" --overwrite
    ;;
  *)
    echo "Unknown MSA tool: $TOOL (expected 'witch' or 'twilight')" >&2
    exit 1
    ;;
esac

rm -rf "$WORKDIR"
