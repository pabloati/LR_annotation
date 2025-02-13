rule busco_run:
    input:
        genome = "data/genome.fasta"
    output:
        directory("busco_output")
    conda:
        "envs/busco.yaml"  
    params:
        busco_dir = "data/busco",
        lineage = config.optional.lineage
    threads:
        config.resources.big.cpu
    shell:
        """
        busco -i {input} -o {output} \
            -l {params.busco_dir} -m genome --augustus \
            -c {threads}
        """

rule busco_gather:
    input:
        "busco_output"  
    output:
        "augustus_model/busco_genes.faa" # Propper file
    params:
        lineage = config.optional.lineage
        gene_type = "single"
    script:
        "scripts/busco_complete_aa.py"

rule clustering_busco_genes:
    input:
        "augustus_model/busco_genes.faa"
    output:
        "augustus_model/cdhit.lst"
    conda:
        "envs/busco.yaml"
    shell:
        """
        cd-hit -o complete_buscos.cdhit -c 0.8 -i {input} -p 1 -d 0 -T 4 -M 48000
        grep ">" complete_buscos.cdhit | cut -f2 -d">" | cut -f1 > {output}
        """

rule concatenate_gff:
    input:
        "augustus_model/chdit.lst"
        "busco_output"
    output:
        directory("augustus_model/busco_genes") # TODO: This has to be a directory 
    conda:
        "envs/busco.yaml"
    params:
        lineage = config.optional.lineage,
        gene_type = "single"
    script:
        "scripts/concatenate_GFF.py"


rule gtf2genbank:
    input:
        genome = config.mandatory.genome,
        gff = "augustus_model/busco_genes/busco_genes.gff"
    output:
        "augustus_model/busco_genes/busco_genes.gb"
    conda:
        "envs/busco.yaml"
    params:
        flanking_reigion = config.optional.flanking_region
    shell:
        """
        gff2gbSmallDNA.pl {input.gff} {input.genome} {params.flanking_region} {output}
        """


rule generate_subsets:
    input:


rule augustus_training:
    input:
