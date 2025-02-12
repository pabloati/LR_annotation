
rule busco_run:
    input:
        genome = "data/genome.fasta"
    output:
        "busco_output"
    conda:
        "envs/busco.yaml"  
    params:
        busco_dir = "data/busco"
    shell:
        """
        busco -i {input} -o {output} \
            -l {params.busco_dir} -m genome
        """

rule busco_to_augustus:
    input:
        busco_output = "busco_output"  
    output:
        "augustus_model/busco_genes.faa" # Propper file
    script:
        "scripts/buco_complete_aa.py"

rule clustering_busco_genes:
    input:
        busco_genes = "augustus_model/busco_genes.faa"
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
    output:
        directory("augustus_model/busco_genes") # TODO: This has to be a directory 
    conda:
        "envs/busco.yaml"
    script:
        "scripts/concatenate_GFF.py"