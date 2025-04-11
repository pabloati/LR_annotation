
if config.isoseq.subreads:
    rule consensus_calling:
        input:
            config.required.input
        output:
            bam = "results/isoannot/input.consensus_calling.bam",
            report = "results/isoannot/input.consensus_calling.report.txt"
        conda:
            f"{dir.envs}/isoseq.yaml"
        shell:
            """
            css {input} {output.bam} --report-file {output.report}
            """
input_type = get_input_type(config.required.input)
if input_type == "bam":
    rule lima:
        input:
            config.required.input
        output:
            name = os.path.join(dir.out.isoseq_lima,"fl.bam"),
            json = os.path.join(dir.out.isoseq_lima,"fl.json")
        conda:
            f"{dir.envs}/isoseq.yaml"
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
            """

    rule lima_renaming:
        input:
            os.path.join(dir.out.isoseq_lima,"fl.json")
        output:
            expand(os.path.join(dir.out.isoseq_lima,"{sample}","{sample}.fl.bam"),sample=samples),
        params:
            samples = config.isoseq.primers_to_samples
        conda:
            f"{dir.envs}/isoseq.yaml"
        threads:
            config.resources.small.cpus,
        resources:
            slurm_extra = f"'--qos={config.resources.small.qos}'",
            cpus_per_task = config.resources.small.cpus,
            mem = config.resources.small.mem,
            runtime =  config.resources.small.time
        log:
            os.path.join(dir.logs,"lima_renaming.log")
        script:
            f"{dir.scripts}/rename_lima_output.py"    
        

    rule refine:
        input:
            lima = os.path.join(dir.out.isoseq_lima,"{sample}","{sample}.fl.bam")
        output:
            os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
        conda:
            f"{dir.envs}/isoseq.yaml"
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
            isoseq refine --require-polya {input.lima} {params.primers} {output} &> {log}
            """

    rule cluster:
        input:
            os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
        output:
            bam=os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam"),
        conda:
            f"{dir.envs}/isoseq.yaml"
        log:
            os.path.join(dir.logs,"isoseq_cluster_{sample}.log")
        threads:
            config.resources.small.cpus,
        resources:
            slurm_extra = f"'--qos={config.resources.medium.qos}'",
            cpus_per_task = config.resources.medium.cpus,
            mem = config.resources.small_bigMem.mem,
            runtime =  config.resources.medium.time
        shell:
            """
            isoseq cluster2 {input} {output.bam} &> {log}
            """

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

elif input_type == "fasta":
    rule IsoQuant:
        input:
            reads = config.required.input
            genome = config.required.genome
        output:
            expand(os.path.join(dir.out.isoseq_isoquant,"{sample}","{sample}.transcript_models.gtf"),sample=samples),
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
            os.path.join(dir.logs,"isoseq_isoquant_{sample}.log")
        shell:
            """
            isoquant.py --reference {input.genome} --fastq {input.reads} \
                --output {dir.out.isoseq_isoquant} --prefix {wildcards.sample} \
                --threads {threads} --data_type {params.data_type} --force 2> {log}
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
        bam = os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam"),
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

rule collapse_isoforms:
    input:
        mapped = os.path.join(dir.out.isoseq_mapping,"{sample}.mapping_pbmm2.bam"),
        flnc = os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
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
        isoseq collapse --do-not-collapse-extra-5exons {input.mapped} {input.flnc} {output} -j {threads} &> {log}
        """
    
# rule filter_transcripts:
#     input:
#         "results/isoannot/{sample}.collapsed"
#     output:
#         "results/isoannot/{sample}.filtered"
#     conda:
#         f"{dir.envs}/isoseq.yaml"
#     resources:
#           slurm_extra = f"'--qos={config.resources.small.qos}'",
#         cpus_per_task = config.resources.small.cpus,
#         mem = config.resources.small.mem,
#         runtime =  config.resources.small.time
#     script:
#         f"{dir.envs}/filter_by_count.py"
