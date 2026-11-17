
import os

#include: "workflows/decompose.smk"
wdr = config["workdir"]
outdir = os.path.join(wdr, "output")
alndir = os.path.join(wdr, "alignments")
rawdir = os.path.join(wdr, "sequences")
bbspec = os.path.join(wdr, "species.txt")
input_bbone = os.path.join(wdr, "backbone.nwk")
bbone = os.path.join(outdir, "backbone.nwk")
trimalndir = os.path.join(outdir, "trimmed")

msa_config = config.get("msa_config", {})
compute_msa = msa_config.get("compute_msa", False)
partition_mode = msa_config.get("partition_mode", "scratch")
backbone_alignment_source = msa_config.get("backbone_alignment", "udance")
userbbalndir = os.path.join(wdr, "backbone_alignments")

if partition_mode not in ("scratch", "extend"):
    raise ValueError("msa_config.partition_mode must be 'scratch' or 'extend', got '%s'" % partition_mode)
if backbone_alignment_source not in ("udance", "user"):
    raise ValueError("msa_config.backbone_alignment must be 'udance' or 'user', got '%s'" % backbone_alignment_source)


def typed_fasta_files(directory, setting):
    # Per-gene files are named <gene>.fna (nucleotide) and/or <gene>.faa (amino acid); when both
    # exist the .faa wins. Returns {gene: (path, "nuc"/"prot")}. chartype still drives every
    # downstream step as a single global setting, so a file of the other type is an error rather
    # than a mistyped alignment fed silently into the rest of the pipeline.
    files = {}
    for f in sorted(os.listdir(directory)):
        base, ext = os.path.splitext(f)
        if ext == ".fna":
            files.setdefault(base, (os.path.join(directory, f), "nuc"))
        elif ext == ".faa":
            files[base] = (os.path.join(directory, f), "prot")
    for path, molecule in files.values():
        if molecule != config["chartype"]:
            raise ValueError(
                "%s: %s is a %s file but config.yaml's chartype is '%s'. "
                "chartype must match the molecule type of your %s input."
                % (setting, path, {"nuc": "nucleotide (.fna)", "prot": "amino-acid (.faa)"}[molecule],
                   config["chartype"], directory)
            )
    return files


USER_BACKBONE_ALNS = {}
if partition_mode == "extend" and backbone_alignment_source == "user":
    if not os.path.isdir(userbbalndir):
        raise ValueError("msa_config.backbone_alignment is 'user' but %s does not exist" % userbbalndir)
    USER_BACKBONE_ALNS = {gene: path for gene, (path, _) in
                          typed_fasta_files(userbbalndir, "msa_config.backbone_alignment").items()}

if compute_msa:
    # RAW_GENE_FILES maps each gene to its chosen raw input path and molecule type; the computed
    # alignment is written as alignments/<gene>.fasta, matching the existing alignments/ convention.
    RAW_GENE_FILES = typed_fasta_files(rawdir, "msa_config.compute_msa")
    GENES = ["%s.fasta" % base for base in RAW_GENE_FILES]

    def compute_msa_input(wildcards):
        base = os.path.splitext(wildcards.gene)[0]
        return RAW_GENE_FILES[base][0]

    def compute_msa_molecule(wildcards):
        base = os.path.splitext(wildcards.gene)[0]
        return RAW_GENE_FILES[base][1]
else:
    GENES = [f for f in os.listdir(alndir) if os.path.isfile(os.path.join(alndir, f))]

TRIMMEDGENES = [os.path.join(outdir, "trimdump", f) for f in GENES]
# Untrimmed per-gene alignments. Partitions take their raw sequences from these (gaps stripped)
# rather than from trimmed/, since trimal -gt drops real residues along with gappy columns.
ALIGNMENTFILES = [os.path.join(alndir, f) for f in GENES]

udance_logpath = os.path.abspath(os.path.join(wdr, "udance.log"))

localrules: all, clean, trimcollect, copybb

rule all:
    input: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])

onstart:
    shell( "if [ -f '{udance_logpath}' ]; then echo '{udance_logpath} already exists."
            "Moving it to {udance_logpath}.bak'; mv {udance_logpath} {udance_logpath}.bak;"
            "else touch {udance_logpath}; fi")

onerror:
    print("Execution failed. Logfile of the execution: ")
    print(udance_logpath)

onsuccess:
    print("Execution is successful. Logfile of the execution: ")
    print(udance_logpath)

rule clean:
    shell:
        """
            rm -r {outdir} 
        """

if compute_msa:
    rule compute_msa:
        input: compute_msa_input
        output: "%s/{gene}" % alndir
        params:
            char=compute_msa_molecule
        resources: cpus=config["resources"]["cores"],
                   mem_mb=config["resources"]["large_memory"]
        benchmark: "%s/benchmarks/compute_msa.{gene}.txt" % outdir
        shell:
            """
            (
            bash uDance/compute_msa.sh {input} {output} {params.char} {resources.cpus} {outdir}/msa_work/{wildcards.gene}
            ) >> {udance_logpath} 2>&1
            """

rule trimtaper:
    input: "%s/{gene}" % alndir
    output: "%s/trimdump/{gene}" % outdir
    params: thr=config["trim_config"]["percent_nongap"]
    benchmark: "%s/benchmarks/trimtaper.{gene}.txt" % outdir
    shell:
        """
            (
            uDance/trimtaper.sh {input} {params.thr} {output}
            ) >> {udance_logpath} 2>&1
        """

rule trimcollect:
    input: TRIMMEDGENES
    output: directory(os.path.join(outdir, "trimmed"))
    shell:
        """
            (
            mv {outdir}/trimdump {output}
            ) >> {udance_logpath} 2>&1
        """

rule mainlines:
    input: trimalndir
    output: os.path.join(outdir, "backbone/0/species.txt")
    params:
            n=config["mainlines_config"]["n"],
            l=config["mainlines_config"]["length"],
            char=config["chartype"],
            bck=config["backbone"]
    resources: mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/mainlines.txt" % outdir
    shell:
        """
            (
            if [ "{params.bck}" == "list" ]; then
                cp {bbspec} {output}
            elif [ "{params.bck}" == "tree" ]; then
                nw_labels -I {input_bbone} > {output}
            elif [ "{params.char}" == "nuc" ]; then  # denovo
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} > {output}
            else
                python run_udance.py mainlines -s {input} -n {params.n} -l {params.l} -p > {output}
            fi
            ) >> {udance_logpath} 2>&1
        """

checkpoint prepbackbonegenes:
    input: s=os.path.join(outdir, "backbone/0/species.txt"), a=trimalndir, raw=ALIGNMENTFILES
    output: touch(os.path.join(outdir,"backbone/0/done.txt"))
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/prepbackbonegenes.txt" % outdir
    params: sub=config["prep_config"]["sublength"],
            frag=config["prep_config"]["fraglength"],
            char=config["chartype"],
            # in tree mode the backbone partition is only needed as the source alignment that
            # extend mode adds each partition onto; its gene trees are never built
            skip="yes" if config["backbone"] == "tree" and not
                 (partition_mode == "extend" and backbone_alignment_source == "udance") else "no"
    shell:
         # python calling shell calling python. looks terrible but
         # this way we are avoiding forking in snakemake main process
        '''
            (
            if [ "{params.skip}" == "yes" ]; then
                touch {outdir}/backbone/0/done.txt
            else
                python -c  "import multiprocessing as mp; \
                            mp.set_start_method('fork'); \
                            from uDance.prep_partition_alignments import prep_partition_alignments; \
                            prep_partition_alignments('{alndir}', \
                                          '{params.char}' == 'prot', \
                                          ['{input.s}'], \
                                          {resources.cpus}, \
                                          {params.sub}, \
                                          {params.frag})"
            fi
            ) >> {udance_logpath} 2>&1
        '''
        # with open(params.logpath, "a") as log_file:
        #     print(params.logpath)



def aggregate_refine_bb_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.prepbackbonegenes.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{j}/seq.fa"))
    return ["%s/backbone/0/%s/bestTree.nwk" % (outdir,j) for j in wc.j]

if config["backbone"] != "tree":
    rule refine_bb:
        input: aggregate_refine_bb_input
        output: bbone = bbone
        params:
            method=config["infer_config"]["method"],
            c=config["refine_config"]["contract"],
            occup = config["refine_config"]["occupancy"],
            ol=config["refine_config"]["outlier_sizelimit"],
            od=config["refine_config"]["outlier_difference"]
        resources: cpus=config["resources"]["cores"],
                   mem_mb=config["resources"]["large_memory"]
        benchmark: "%s/benchmarks/refine_copy_bb.txt" % outdir
        shell:
            '''
                (
                python run_udance.py refine -p {outdir}/backbone/0 -m {params.method} -M {resources.mem_mb} -c {params.c} -o {params.occup} -T {resources.cpus} -l {params.ol} -d {params.od}
                nw_reroot -d {outdir}/backbone/0/astral_output.incremental.nwk > {output}
                ) >> {udance_logpath} 2>&1
            '''
else:
    rule copybb:
        input: os.path.join(outdir,"backbone/0/done.txt")
        output: bbone = bbone
        benchmark: "%s/benchmarks/refine_copy_bb.txt" % outdir
        shell:
            '''
                cp {input_bbone} {output.bbone}
            '''

rule placement_prep:
    input: b = bbone,
           ind = trimalndir
    output: aln = os.path.join(outdir,"placement/backbone.fa"),
            qry = os.path.join(outdir,"placement/query.fa"),
            tre = os.path.join(outdir,"placement/backbone.tree")
    params: char=config["chartype"],
            filtering=config["backbone_filtering"],
            f=config["apples_config"]["filter"],
            m=config["apples_config"]["method"],
            b=config["apples_config"]["base"],
            v=config["apples_config"]["overlap"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/placement_prep.txt" % outdir
    shell:
        """
            (
            bash uDance/create_concat_alignment.sh {input.ind} {input.b} {outdir} {params.char} {resources.cpus} {params.f} {params.m} {params.b} {params.v} {params.filtering}
            ) >> {udance_logpath} 2>&1            
        """

rule placement:
    input: aln = os.path.join(outdir,"placement/backbone.fa"),
           qry = os.path.join(outdir,"placement/query.fa"),
           tre = os.path.join(outdir,"placement/backbone.tree")
    output: j=os.path.join(outdir,"placement.jplace")
    params: f=config["apples_config"]["filter"],
            m=config["apples_config"]["method"],
            b=config["apples_config"]["base"],
            v=config["apples_config"]["overlap"],
            char=config["chartype"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/placement.txt" % outdir
    log: out=os.path.join(outdir,"placement/apples2.out"), err=os.path.join(outdir,"placement/apples2.err")
    shell:
        """
            (
            export MKL_NUM_THREADS=1
            export NUMEXPR_NUM_THREADS=1
            export OMP_NUM_THREADS=1
            if [ "{params.char}" == "nuc" ]; then
                run_apples.py --exclude -s {input.aln} -q {input.qry} -T {resources.cpus} -V {params.v} \
                -t {input.tre} -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
            else
                run_apples.py --exclude -p -s {input.aln} -q {input.qry} -T {resources.cpus} -V {params.v} \
                 -t {input.tre} -f {params.f} -m {params.m} -b {params.b} -o {output.j} > {log.out} 2> {log.err}
            fi
            ) >> {udance_logpath} 2>&1            
        """


checkpoint decompose:
    input: j=os.path.join(outdir,"placement.jplace"), ind=trimalndir, raw=ALIGNMENTFILES
    output: cst=os.path.join(outdir,"udance/color_spanning_tree.nwk")
    params:
        size=config["prep_config"]["cluster_size"],
        method=config["infer_config"]["method"],
        edg=config["prep_config"]["edge_thr"],
        sub=config["prep_config"]["sublength"],
        frag=config["prep_config"]["fraglength"],
        pra=config["prep_config"]["pruneafter"],
        mps=config["prep_config"]["min_placements"],
        char=config["chartype"]

    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/decompose.txt" % outdir
    shell:
        """
            (
            cp {input.j} {outdir}/udance
            if [ "{params.size}" == "auto" ]; then
                numq=`grep ">" {outdir}/placement/query.fa  | wc -l`
                numbb=`grep ">" {outdir}/placement/backbone.fa  | wc -l`
                clustsz=`python -c "import math; print(int(min(2500,round(8*math.sqrt($numbb)+3*math.sqrt($numq),-2))))"`
                echo "Cluster size is set to $clustsz (automatically)"
            else
                clustsz="{params.size}"
                echo "Cluster size is set to $clustsz (user-choice)"
            fi
            if [ "{params.char}" == "nuc" ]; then
                python run_udance.py decompose -s {input.ind} --sequence-dir {alndir} -o {outdir}/udance -t $clustsz -j {input.j} \
                -m {params.method} -T {resources.cpus} -l {params.sub} -f {params.frag} -e {params.edg} \
                --minplacements {params.mps}
            else
                python run_udance.py decompose -p -s {input.ind} --sequence-dir {alndir} -o {outdir}/udance -t $clustsz -j {input.j} \
                -m {params.method} -T {resources.cpus} -l {params.sub} -f {params.frag} -e {params.edg} \
                --minplacements {params.mps}
            fi
            python prune_similar.py -T {resources.cpus} -o {outdir}/udance -S {params.pra}
            if [  -f {outdir}/udance/dedupe_map.txt ]; then 
                cat {outdir}/udance/dedupe_map.txt > {outdir}/dedupe_map.txt 
            fi 
            if [ -f {outdir}/rm_map.txt ]; then 
                cat {outdir}/rm_map.txt >> {outdir}/dedupe_map.txt 
            fi
            ) >> {udance_logpath} 2>&1
        """

def partition_msa_backbone(wildcards):
    # The alignment extend mode adds this partition's queries onto, or [] to align from scratch.
    if partition_mode != "extend":
        return []
    if backbone_alignment_source == "user":
        return [USER_BACKBONE_ALNS[wildcards.gene]] if wildcards.gene in USER_BACKBONE_ALNS else []
    if wildcards.stage == "backbone":
        # this partition is the source alignment itself
        return []
    checkpoints.prepbackbonegenes.get()
    if os.path.exists(os.path.join(outdir, "backbone/0", wildcards.gene, "seq.fa")):
        return [os.path.join(outdir, "backbone/0", wildcards.gene, "aln.fa")]
    return []

# per-partition MSA: aligns each partition's raw sequences of one gene, then trims them
rule partition_msa:
    input:
        seq="%s/{stage}/{cluster}/{gene}/seq.fa" % outdir,
        backbone=partition_msa_backbone
    output:
        "%s/{stage}/{cluster}/{gene}/aln.fa" % outdir
    params:
        mode=partition_mode,
        bb=lambda wildcards: (partition_msa_backbone(wildcards) or ["none"])[0],
        char=config["chartype"],
        trimthr=config["trim_config"]["percent_nongap"],
        work=lambda wildcards: os.path.join(outdir, "msa_work", wildcards.stage, wildcards.cluster, wildcards.gene)
    log: "%s/{stage}/{cluster}/{gene}/msa.log" % outdir
    benchmark: "%s/{stage}/{cluster}/{gene}/msa_benchmark.txt" % outdir
    shell:
        '''
            bash uDance/partition_msa.sh {input.seq} {output} {params.mode} {params.bb} \
                {params.char} {params.trimthr} {threads} {params.work} > {log} 2>&1
        '''

# phy inf
rule genetreeinfer:
    input:
        "%s/{stage}/{cluster}/{gene}/aln.fa" % outdir
    output:
        "%s/{stage}/{cluster}/{gene}/bestTree.nwk" % outdir
    params:
          c=config["chartype"],
          s=config["infer_config"]["numstart"],
          thrd=config["infer_config"]["numthread"],
          t=config["infer_config"]["method"]
    benchmark: "%s/{stage}/{cluster}/{gene}/benchmark.txt" % outdir
    shell:
        '''
            # many instances of this rule may run simultaneously. To reduce the IO overhead, we "sponge" the output
            # before appending to udance_logpath
            source uDance/mysponge.sh
            (
            bash uDance/process_a_marker.sh {input} {params.c} {params.s} {params.t} {params.thrd}
            ) 2>&1 | mysponge #>> {udance_logpath} 
        '''



def aggregate_refine_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "%s/{j}/seq.fa" % wildcards.cluster))
    return ["%s/udance/%s/%s/bestTree.nwk" % (outdir, wildcards.cluster,j) for j in wc.j]


rule refine:
    input: aggregate_refine_input
    output: expand("%s/udance/{{cluster}}/astral_output.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: method=config["infer_config"]["method"],
            c=config["refine_config"]["contract"],
            occup=config["refine_config"]["occupancy"],
            ol=config["refine_config"]["outlier_sizelimit"],
            od=config["refine_config"]["outlier_difference"]
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/refine.{cluster}.txt" % outdir

    shell:
        """
            (
            python run_udance.py refine -p {outdir}/udance/{wildcards.cluster} -m {params.method} -M {resources.mem_mb} -c {params.c} -o {params.occup} -T {resources.cpus} -l {params.ol} -d {params.od}
            ) >> {udance_logpath} 2>&1
        """

rule blinference:
    input: expand("%s/udance/{{cluster}}/astral_output.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    output: expand("%s/udance/{{cluster}}/astral_output.{approach}.nwk.bl" % outdir, approach=["incremental", "updates"])
    resources: cpus=config["resources"]["cores"],
               mem_mb=config["resources"]["large_memory"]
    benchmark: "%s/benchmarks/blinference.{cluster}.txt" % outdir
    shell:
        '''
            pwdd=`pwd`
            for approach in incremental updates; do
                if [ -f  {outdir}/udance/{wildcards.cluster}/skip_partition ] ; then
                    cp {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl
                else
                    java -Xmx{resources.mem_mb}M -Djava.library.path=$pwdd/uDance/tools/ASTRAL/lib/ -jar $pwdd/uDance/tools/ASTRAL/astralmp.5.17.2.jar \
                        -q {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk \
                        -i {outdir}/udance/{wildcards.cluster}/astral_input.trees \
                        -o {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl \
                        -C -T {resources.cpus} -u > {outdir}/udance/{wildcards.cluster}/astral.$approach.log.bl 2>&1
                    mv {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl.falsesupport
                    python uDance/transfer_supports.py {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl.falsesupport > {outdir}/udance/{wildcards.cluster}/astral_output.$approach.nwk.bl
                fi
            done
        '''

def aggregate_stitch_input(wildcards):
    checkpoint_output = os.path.dirname(checkpoints.decompose.get(**wildcards).output[0])
    wc = glob_wildcards(os.path.join(checkpoint_output, "{i}/species.txt"))
    if config["refine_config"]["infer_branchlen"] in [False, "False"]:
        return [f"%s/udance/%s/astral_output.%s.nwk" % (outdir, i, j) for i in wc.i for j in ["incremental", "updates"]]
    else:
        return [f"%s/udance/%s/astral_output.%s.nwk.bl" % (outdir, i, j) for i in wc.i for j in ["incremental", "updates"]]


rule stitch:
    input: aggregate_stitch_input
    output: expand("%s/udance.{approach}.nwk" % outdir, approach=["incremental", "updates"])
    params: b = config["refine_config"]["infer_branchlen"]
    benchmark: "%s/benchmarks/stitch.txt" % outdir
    shell:
        """
           (
            if [[ "{params.b}" == "False" ]] ; then
                python run_udance.py stitch -o {outdir}/udance
            else
                python run_udance.py stitch -o {outdir}/udance -b
            fi
            cp {outdir}/udance/udance.*.nwk {outdir}
           ) >> {udance_logpath} 2>&1
        """
