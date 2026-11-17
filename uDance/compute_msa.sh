#!/usr/bin/env bash
set -euo pipefail

INP=$1      # raw unaligned fasta
OUT=$2      # aligned output
CHAR=$3     # nuc | prot (uDance's existing chartype values)
THR=$4      # cpu threads
WORKDIR=$5  # scratch dir, unique per call

mkdir -p "$WORKDIR"

molecule="dna"
[ "$CHAR" == "prot" ] && molecule="amino"
witch-msa -i "$INP" -d "$WORKDIR" -o aligned.fasta -t "$THR" --molecule "$molecule"
mv "$WORKDIR/aligned.fasta" "$OUT"

rm -rf "$WORKDIR"
