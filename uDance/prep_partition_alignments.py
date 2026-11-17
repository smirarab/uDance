import multiprocessing as mp

from uDance.PoolAlignmentWorker import PoolAlignmentWorker
from uDance.fasta2dic import readfq

from os import listdir
from os.path import isfile, join, splitext


# seqdir may hold aligned or unaligned FASTA files: gaps are stripped either way, because each
# partition is re-aligned on its own by the partition_msa rule.
def prep_partition_alignments(
    seqdir, protein_flag, species_path_list, num_thread, subalignment_length, fragment_length
):
    if protein_flag:
        valid_translation = str.maketrans('BJOUXZ', '-' * 6)
    else:
        valid_translation = str.maketrans('BDEFHIJKLMNOPQRSUVWXYZ', '-' * 22)
    only_files = [f for f in listdir(seqdir) if isfile(join(seqdir, f)) and not f.startswith('.')]
    all_scripts = []
    for aln in only_files:
        seq_input_file = join(seqdir, aln)
        basename = splitext(aln)[0]
        with open(seq_input_file) as f:
            seq_dict = {name: seq.upper().replace('-', '').replace('.', '') for name, seq, _ in readfq(f)}
        alignment_worker = PoolAlignmentWorker()
        alignment_worker.set_class_attributes(
            subalignment_length, fragment_length, seq_dict, valid_translation, basename
        )
        pool = mp.Pool(num_thread)
        scripts = pool.map(alignment_worker.worker, species_path_list)
        pool.close()
        pool.join()
        valid_scripts = [s for s in scripts if s is not None]
        all_scripts += valid_scripts
    return all_scripts
