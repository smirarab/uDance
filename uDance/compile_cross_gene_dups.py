#!/usr/bin/env python
"""Phase 2 of cross-gene deduplication.

Reads per-gene seqkit duplicate lists, compiles cross-gene fingerprints,
identifies genomes that are identical across ALL genes they appear in
(same gene set required), and writes per-gene dupmap/remove files.

Usage:
    python -I compile_cross_gene_dups.py <dup_dir> <output_dir> <input_aln_dir> <backbone_taxa_file>

Where:
    dup_dir:            directory with <gene>.dup.txt files from seqkit rmdup -D
    output_dir:         directory to write dupmap.txt, remove.txt, and dedup_map.txt
    input_aln_dir:      original alignment directory (to get full genome lists per gene)
    backbone_taxa_file: file with backbone taxon names (one per line)
"""

import os
import sys
from collections import defaultdict


def parse_seqkit_dupfile(filepath):
    """Parse seqkit rmdup -D output.

    Format: count<TAB>name1, name2, name3
    Returns list of lists, each a group of duplicate genome names.
    """
    groups = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t', 1)
            if len(parts) < 2:
                continue
            names = [n.strip() for n in parts[1].split(', ')]
            if len(names) >= 2:
                groups.append(names)
    return groups


def get_genome_names(fasta_path):
    """Extract genome names from fasta file (fast, no sequence reading)."""
    names = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                names.append(line[1:].split()[0].strip())
    return names


def main():
    if len(sys.argv) != 5:
        print('Usage: compile_cross_gene_dups.py <dup_dir> <output_dir> <input_aln_dir> <backbone_taxa_file>',
              file=sys.stderr)
        sys.exit(1)

    dup_dir = sys.argv[1]
    output_dir = sys.argv[2]
    input_aln_dir = sys.argv[3]
    backbone_taxa_file = sys.argv[4]

    os.makedirs(output_dir, exist_ok=True)

    with open(backbone_taxa_file) as f:
        backbone_taxa = set(line.strip() for line in f if line.strip())

    # Phase 2a: build per-genome signatures from seqkit dup files and full genome lists
    # genome_sigs[genome][gene] = representative name (first in seqkit dup group, or self if unique)
    genome_sigs = defaultdict(dict)

    gene_names = []
    for fname in sorted(os.listdir(input_aln_dir)):
        fpath = os.path.join(input_aln_dir, fname)
        if not os.path.isfile(fpath):
            continue
        gene_names.append(fname)

        all_genomes = get_genome_names(fpath)

        # parse seqkit dup file if it exists
        dup_file = os.path.join(dup_dir, fname + '.dup.txt')
        dup_groups = []
        genome_to_rep = {}
        if os.path.isfile(dup_file):
            dup_groups = parse_seqkit_dupfile(dup_file)
            for group in dup_groups:
                rep = group[0]
                for name in group:
                    genome_to_rep[name] = rep

        for genome in all_genomes:
            rep = genome_to_rep.get(genome, genome)
            genome_sigs[genome][fname] = rep

    print('Read %d genomes across %d genes.' % (len(genome_sigs), len(gene_names)))

    # Phase 2b: group by fingerprint
    fingerprint_groups = defaultdict(list)
    for genome, sigs in genome_sigs.items():
        fp = tuple(sorted(sigs.items()))
        fingerprint_groups[fp].append(genome)

    # Phase 2c: identify duplicate groups (2+ members with identical fingerprint)
    dup_groups = []
    for members in fingerprint_groups.values():
        if len(members) >= 2:
            # representative selection: prefer backbone taxa
            backbone_members = sorted(m for m in members if m in backbone_taxa)
            query_members = sorted(m for m in members if m not in backbone_taxa)
            if backbone_members:
                ordered = backbone_members + query_members
            else:
                ordered = sorted(members)
            dup_groups.append(ordered)

    total_dups = sum(len(g) - 1 for g in dup_groups)
    print('Found %d duplicate groups (%d genomes to remove).' % (len(dup_groups), total_dups))

    if not dup_groups:
        print('No cross-gene duplicates found.')
        return

    # collect all duplicates (non-representatives) to remove
    all_remove = set()
    for group in dup_groups:
        for name in group[1:]:
            all_remove.add(name)

    # Phase 2d: write outputs
    # global dedup_map.txt
    dedup_map_path = os.path.join(output_dir, 'dedup_map.txt')
    with open(dedup_map_path, 'w') as f:
        for group in sorted(dup_groups):
            f.write('\t'.join(group) + '\n')

    # per-gene dupmap.txt and remove.txt
    for gene in gene_names:
        gene_genomes = set(get_genome_names(os.path.join(input_aln_dir, gene)))

        gene_dup_groups = []
        gene_remove = []
        for group in dup_groups:
            members_here = [m for m in group if m in gene_genomes]
            if len(members_here) >= 2:
                gene_dup_groups.append(members_here)
                for m in members_here[1:]:
                    gene_remove.append(m)

        if gene_dup_groups:
            dupmap_path = os.path.join(output_dir, gene + '.dupmap.txt')
            with open(dupmap_path, 'w') as f:
                for group in gene_dup_groups:
                    f.write('\t'.join(group) + '\n')

        if gene_remove:
            remove_path = os.path.join(output_dir, gene + '.remove.txt')
            with open(remove_path, 'w') as f:
                for name in gene_remove:
                    f.write(name + '\n')

    print('Wrote dedup_map.txt and per-gene dupmap/remove files to %s' % output_dir)


if __name__ == '__main__':
    main()
