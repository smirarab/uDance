#!/usr/bin/env bash
set -euo pipefail

SEQ=$1       # partition's raw sequences for one gene (<stage>/<cluster>/<gene>/seq.fa)
OUT=$2       # aligned + trimmed output (<stage>/<cluster>/<gene>/aln.fa), read by genetreeinfer
MODE=$3      # scratch | extend
BACKBONE=$4  # backbone alignment for extend mode, or "none"
CHAR=$5      # nuc | prot
TRIMTHR=$6   # trim_config.percent_nongap
THR=$7       # cpu threads
WORKDIR=$8   # scratch dir, unique per partition x gene

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"

align_from_scratch() {
    bash uDance/compute_msa.sh "$SEQ" "$WORKDIR/msa.fa" "$CHAR" "$THR" "$WORKDIR/scratch"
}

if [ "$MODE" == "extend" ] && [ "$BACKBONE" != "none" ]; then
    read -r NB NQ < <(python -m uDance.split_partition "$SEQ" "$BACKBONE" "$WORKDIR/backbone.fa" "$WORKDIR/queries.fa")
    if [ "$NB" -lt 4 ]; then
        # too few backbone rows to build HMMs from
        echo "partition_msa: only $NB backbone rows in $SEQ, aligning it from scratch instead"
        align_from_scratch
    elif [ "$NQ" -eq 0 ]; then
        cp "$WORKDIR/backbone.fa" "$WORKDIR/msa.fa"
    else
        molecule="dna"
        [ "$CHAR" == "prot" ] && molecule="amino"
        witch-msa -b "$WORKDIR/backbone.fa" -q "$WORKDIR/queries.fa" -d "$WORKDIR/witch" \
            -o aligned.fasta -t "$THR" --molecule "$molecule"
        mv "$WORKDIR/witch/aligned.fasta" "$WORKDIR/msa.fa"
    fi
else
    align_from_scratch
fi

# insertion columns may come back lower-case; everything downstream expects upper case
awk '/^>/{print; next} {print toupper($0)}' "$WORKDIR/msa.fa" > "$WORKDIR/msa_upper.fa"
bash uDance/trimtaper.sh "$WORKDIR/msa_upper.fa" "$TRIMTHR" "$OUT"

rm -rf "$WORKDIR"
