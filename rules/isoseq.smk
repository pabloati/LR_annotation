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
    shell:
        """
        lima {input} {params.primers} {output.name} \
            --isoseq --peek-guess  --split-subdirs --overwrite-biosample-names \
            --split-named --biosample-csv {params.samples}
        touch {output}
        """

rule lima_renaming:
    input:
        os.path.join(dir.out.isoseq_lima,"fl.json")
    output:
        touch(os.path.join(dir.out.isoseq_lima,"renamed.done")),
        expand(os.path.join(dir.out.isoseq_lima,"{sample}","{sample}.fl.bam"),sample=samples)
    params:
        samples = config.isoseq.primers_to_samples
    conda:
        f"{dir.env}/isoseq.yaml"
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
    params:
        primers = config.isoseq.primers
    shell:
        """
        isoseq3 refine --require-polya {input.lima} {params.primers} {output}
        """

rule cluster:
    input:
        os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
    output:
        bam=os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam"),
    conda:
        f"{dir.env}/isoseq.yaml"
    shell:
        """
        isoseq3 cluster {input} {output.bam} --use-qvs
        """

rule bam2fastq:
    input:
        os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam")
    output:
        os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.fastq")
    conda:
        f"{dir.env}/isoseq.yaml"
    shell:
        """
        samtools fastq {input} > {output}
        """

rule mapping_reads:
    input:
        reads=os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.fastq"),
        genome = config.required.genome
    output:
        os.path.join(dir.out.isoseq_mapping,"{sample}.mapping.sam"),
    conda:
        f"{dir.env}/minimap2.yaml"
    shell:
        """
        minimap2 -ax splice  -uf --secondary=no -C5 {input.genome} {input.reads} | sort -k 3,3 -k 4,4n > {output}
        """

rule collapse_isoforms:
    input:
        os.path.join(dir.out.isoseq_mapping,"{sample}.mapping.sam"),
        isoforms=os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam"),
    output:
        "results/isoannot/{sample}.collapsed"
    conda:
        f"{dir.env}/isoseq.yaml"
    script:
        "scripts/collapse_isoforms_by_sam.py"
    
rule count_reads:
    input:
        "results/isoannot/{sample}.collapsed"
    output:
        "results/isoannot/polished.cluster_report.csv"
    conda:
        f"{dir.env}/isoseq.yaml"
    script:
        "scripts/get_abundance_post_collapse.py"

rule filter_transcripts:
    input:
        "results/isoannot/{sample}.collapsed"
    output:
        "results/isoannot/{sample}.filtered"
    conda:
        f"{dir.env}/isoseq.yaml"
    script:
        "scripts/filter_by_count.py"