rule consensus_calling:
    input:
        "data/consensus_calling/{sample}.bam"
    output:
        bam = "results/isoannot/{sample}.consensus_calling.bam",
        report = "results/isoannot/{sample}.consensus_calling.report.txt"
    conda:
        "envs/isoseq.yaml"
    shell:
        """
        css {input} {output.bam} --report-file {output.report}
        """

rule lima:
    input:
        "results/isoannot/{sample}.consensus_calling.bam"
    output:
        "results/isoannot/{sample}.fl.bam"
    conda:
        "envs/isoseq.yaml"
    params:
        primers = "data/primers.fasta"
    shell:
        """
        lima --dump-clips --isoseq {input} {params.primers} {output}
        """

rule refine:
    input:
        "results/isoannot/{sample}.fl.bam"
    output:
        "results/isoannot/{sample}.refine.bam"
    conda:
        "envs/isoseq.yaml"
    params:
        primers = "data/primers.fasta"
    shell:
        """
        isoseq3 refine --require-polyA {input} {params.primers} {output}
        """

rule cluster:
    input:
        "results/isoannot/{sample}.refine.bam"
    output:
        bam="results/isoannot/{sample}.cluster.bam",
        reads="results/isoannot/{sample}.cluster.fastq"
    conda:
        "envs/isoseq.yaml"
    shell:
        """
        isoseq3 cluster {input} {output.bam} --verbose --use-qvs
        """

rule mapping_reads:
    input:
        reads="data/reads/{sample}.cluster.fastq",
        genome="data/genome.fasta"
    output:
        "results/isoannot/{sample}.mapping.sam"
    conda:
        "envs/isoseq.yaml"
    shell:
        """
        minimap2 -ax splice  -uf --secondary=no -C5 {input.genome} {input.reads} | sort -k 3,3 -k 4,4n > {output}
        """

rule collapse_isoforms:
    input:
        map_reads="results/isoannot/{sample}.mapping.sam",
        isoforms="results/isoannot/{sample}.cluster.bam"
    output:
        "results/isoannot/{sample}.collapsed"
    conda:
        "envs/isoseq.yaml"
    script:
        "scripts/collapse_isoforms_by_sam.py"
    
rule count_reads:
    input:
        "results/isoannot/{sample}.collapsed"
    output:
        "results/isoannot/polished.cluster_report.csv"
    conda:
        "envs/isoseq.yaml"
    script:
        "scripts/get_abundance_post_collapse.py"

rule filter_transcripts:
    input:
        "results/isoannot/{sample}.collapsed"
    output:
        "results/isoannot/{sample}.filtered"
    conda:
        "envs/isoseq.yaml"
    script:
        "scripts/filter_by_count.py"