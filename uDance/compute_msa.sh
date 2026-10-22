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
    # twilight lives in its own conda env (udance_twilight), not the main "udance" env:
    # its bioconda build needs libzlib>=1.3.1, which conflicts with gappa=0.7.1's
    # zlib<1.3.0a0 pin in the main env (see install.sh for the full explanation).
    TWILIGHT_BIN="$(conda info --base)/envs/udance_twilight/bin/twilight"
    if [ ! -x "$TWILIGHT_BIN" ]; then
        echo "ERROR: $TWILIGHT_BIN not found or not executable." >&2
        echo "Run install.sh to create the udance_twilight conda environment first." >&2
        exit 1
    fi
    # Step 1: quick guide tree via MAFFT's --treeout side-output.
    mafft --retree 0 --treeout --reorder --thread "$THR" "$INP" > "$WORKDIR/mafft.out"
    python3 uDance/mafft2nwk.py "$INP.tree" "$INP" "$WORKDIR/guide.nwk"
    rm -f "$INP.tree" "$WORKDIR/mafft.out"
    # Step 2: TWILIGHT does the actual progressive alignment along that guide tree.
    "$TWILIGHT_BIN" -i "$INP" -t "$WORKDIR/guide.nwk" -o "$OUT" -C "$THR" --type "$type" --overwrite
    ;;
  *)
    echo "Unknown MSA tool: $TOOL (expected 'witch' or 'twilight')" >&2
    exit 1
    ;;
esac

rm -rf "$WORKDIR"
