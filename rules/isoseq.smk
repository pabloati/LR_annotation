with open(config.isoseq.primers_to_samples) as f:
    samples = [line.strip().split(",")[1] for line in f.readlines()[1:]]
    barcodes = [line.strip().split(",")[0] for line in f.readlines()[1:]]
if config.isoseq.subreads:
    rule consensus_calling:
        input:
            config.required.input
        output:
            bam = "results/isoannot/input.consensus_calling.bam",
            report = "results/isoannot/input.consensus_calling.report.txt"
        conda:
            f"{dir.env}/isoseq.yaml"
        shell:
            """
            css {input} {output.bam} --report-file {output.report}
            """

rule lima:
    input:
        config.required.input
    output:
        name = os.path.join(dir.out.isoseq_lima,"fl.bam"),
        json = os.path.join(dir.out.isoseq_lima,"fl.json")
    conda:
        f"{dir.env}/isoseq.yaml"
    params:
        primers = config.isoseq.primers,
        samples = config.isoseq.primers_to_samples
    log:
        os.path.join(dir.logs,"lima_demultiplexing.log")
    threads:
        config.resources.big.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.big.qos}'",
        cpus_per_task = config.resources.big.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    shell:
        """
        lima {input} {params.primers} {output.name} \
            --isoseq --peek-guess  --split-subdirs --overwrite-biosample-names \
            --split-named --biosample-csv {params.samples} &> {log}
        touch {output}
        """

rule lima_renaming:
    input:
        os.path.join(dir.out.isoseq_lima,"fl.json")
    output:
        touch(os.path.join(dir.out.isoseq_lima,"renamed.done")),
        expand(os.path.join(dir.out.isoseq_lima,"{sample}","{sample}.fl.bam"),sample=samples),
    params:
        samples = config.isoseq.primers_to_samples
    conda:
        f"{dir.env}/isoseq.yaml"
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    script:
        f"{dir.scripts}/rename_lima_output.py"    
    

rule refine:
    input:
        os.path.join(dir.out.isoseq_lima,"renamed.done"),
        lima = os.path.join(dir.out.isoseq_lima,"{sample}","{sample}.fl.bam")
    output:
        os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
    conda:
        f"{dir.env}/isoseq.yaml"
    log:
        os.path.join(dir.logs,"isoseq_refine_{sample}.log")
    params:
        primers = config.isoseq.primers
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        isoseq3 refine --require-polya {input.lima} {params.primers} {output} &> {log}
        """

rule cluster:
    input:
        os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
    output:
        bam=os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam"),
    conda:
        f"{dir.env}/isoseq.yaml"
    log:
        os.path.join(dir.logs,"isoseq_cluster_{sample}.log")
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        isoseq3 cluster {input} {output.bam} --use-qvs &> {log}
        """

rule bam2fastq:
    input:
        os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam")
    output:
        os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.fastq")
    conda:
        f"{dir.env}/minimap2.yaml"
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        samtools fastq {input} > {output}
        """

rule mapping_reads:
    input:
        reads=os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.fastq"),
        genome = config.required.genome
    output:
        os.path.join(dir.out.isoseq_mapping,"{sample}.mapping.bam"),
    conda:
        f"{dir.env}/minimap2.yaml"
    threads:
        config.resources.big.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.big.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    shell:
        """
        minimap2 -ax splice -t {threads} -uf --secondary=no -C5 {input.genome} {input.reads} | samtools sort |  samtools view -bS > {output}
        """

rule collapse_isoforms:
    input:
        mapped = os.path.join(dir.out.isoseq_mapping,"{sample}.mapping.bam"),
        flnc = os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
    output:
        os.path.join(dir.out.isoseq_collapsed,"{sample}","{sample}.collapsed.gff"),
    conda:
        f"{dir.env}/isoseq.yaml"
    log:
        os.path.join(dir.logs,"isoseq_collapse_{sample}.log")
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        isoseq collapse --do-not-collapse-extra-5exons {input.mapped} {input.flnc} {output} -j {threads} &> {log}
        """
    
# rule filter_transcripts:
#     input:
#         "results/isoannot/{sample}.collapsed"
#     output:
#         "results/isoannot/{sample}.filtered"
#     conda:
#         f"{dir.env}/isoseq.yaml"
#     resources:
#           slurm_extra = f"'--qos={config.resources.small.qos}'",
#         cpus_per_task = config.resources.small.cpus,
#         mem = config.resources.small.mem,
#         runtime =  config.resources.small.time
#     script:
#         f"{dir.env}/filter_by_count.py"