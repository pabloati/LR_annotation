


# TODO: Run this rule only if the input is a BAM file and IsoQANT is not used
rule bam2fastq:
    input:
        os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam")
    output:
        os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.fastq")
    conda:
        f"{dir.envs}/minimap2.yaml"
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
        f"{dir.envs}/minimap2.yaml"
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

rule index_genome:
    input:
        genome = config.required.genome
    output:
        os.path.join(dir.tools_index,genome_name,"index.mmi")
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads:
        config.resources.medium.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.medium.qos}'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.medium.time
    shell:
        """
        mkdir -p {dir.tools_index}
        pbmm2 index {input.genome} {output}
        """

rule mapping_reads_pbmm2:
    input:
        bam = config.required.input,
        index = os.path.join(dir.tools_index,genome_name,"index.mmi")
    output:
        os.path.join(dir.out.isoseq_mapping,"{sample}.mapping_pbmm2.bam"),
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads:
        config.resources.big.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.big.qos}'",
        cpus_per_task = config.resources.big.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    log:
        os.path.join(dir.logs,"isoseq_mapping_{sample}.log")
    shell:
        """
        pbmm2 align --preset ISOSEQ --sort {input.index} {input.bam}  {output} &> {log}
        """

# TODO: see how to use the FLNC BAM if we decide to use IsoSeq3
rule collapse_isoforms:
    input:
        mapped = os.path.join(dir.out.isoseq_mapping,"{sample}.mapping_pbmm2.bam"),
        #flnc = os.path.join(config.required.flnc_dir,"{sample}","{sample}.flnc.bam")
    output:
        os.path.join(dir.out.isoseq_collapsed,"{sample}","{sample}.collapsed.gff"),
    conda:
        f"{dir.envs}/isoseq.yaml"
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
        isoseq collapse --do-not-collapse-extra-5exons {input.mapped} {output} -j {threads} &> {log}
        """
