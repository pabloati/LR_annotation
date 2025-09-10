


# TODO: Run this rule only if the input is a BAM file and IsoQANT is not used

rule bam2fastq:
    input:
        config.required.input
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

rule IsoQuant:
    input:
        reads = os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.fastq"),
        genome = config.required.genome
    output:
        os.path.join(dir.out.isoquant,"{sample}","{sample}.transcript_models.gtf"),
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    conda:
        f"{dir.envs}/isoquant.yaml"
    params:
        data_type = config.isoquant.data_type,
    threads:
        config.resources.small.cpus,
    log:
        os.path.join(dir.logs,"isoseq_isoquant.log")
    shell:
        """
        prefix=$(basename {output} .transcript_models.gtf)
        isoquant.py --reference {input.genome} --fastq {input.reads} \
            --output {dir.out.isoquant} --prefix $prefix \
            --threads {threads} --data_type {params.data_type} --force &> {log}
        """