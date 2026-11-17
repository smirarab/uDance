import sys

from uDance.fasta2dic import readfq

# usage: python -m uDance.split_partition <seq.fa> <backbone alignment> <out backbone.fa> <out queries.fa>
# Taxa of the partition that appear in the backbone alignment keep their aligned rows (raw rows
# for them in seq.fa are ignored); every other taxon becomes a query to be added onto them.
# Prints "<number of backbone rows> <number of queries>".

seq_fp, backbone_fp, out_backbone_fp, out_queries_fp = sys.argv[1:5]

with open(seq_fp) as f:
    partition = {name: seq for name, seq, _ in readfq(f)}
with open(backbone_fp) as f:
    backbone = {name: seq.upper() for name, seq, _ in readfq(f) if name in partition}

keep = [i for i in range(len(next(iter(backbone.values()), ''))) if any(s[i] != '-' for s in backbone.values())]
with open(out_backbone_fp, 'w') as f:
    for name, seq in backbone.items():
        f.write('>%s\n%s\n' % (name, ''.join(seq[i] for i in keep)))

queries = {name: seq for name, seq in partition.items() if name not in backbone}
with open(out_queries_fp, 'w') as f:
    for name, seq in queries.items():
        f.write('>%s\n%s\n' % (name, seq))

print(len(backbone), len(queries))
