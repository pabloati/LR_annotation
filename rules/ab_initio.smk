# Snakefile for ab initio gene prediction
#TODO: move to a module
def calculate_subset_size(file_path,test_size):
    count = 0
    try:
        with open(file_path, 'r') as file:
            for line in file:
                count += line.upper().count("LOCUS")
        return count
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return None
    except IOError:
        print(f"Error: Unable to read the file '{file_path}'.")
        return None
    # Subset based on size
    testing_transcripts = count - test_size
    # If the number of transcripts is less than 1000
    if testing_transcripts < 1000:
        subset = [testing_transcripts]
    
    # if the testing_transcripts is between 1000 and 1500
    elif testing_transcripts > 1000 and testing_transcripts < 1500:
        subset = [200, 500, 800, 1000]
    
    # if the testing_transcripts is between 1500 and 2000
    elif testing_transcripts > 1500 and testing_transcripts < 2000:
        subset = [200, 500, 800, 1000, 1500]
    
    # From 2000 increase the number by 1000 to 8000
    elif testing_transcripts > 2000:
        subset = [200, 500, 800, 1000, 1500, 2000]
        while (subset[-1] + 1000) < testing_transcripts and (subset[-1] + 1000) <= 8000:
            subset.append(subset[-1] + 1000)
    return ','.join(map(str, subset))

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
        config.resources.big.cpus
    log:
        "logs/busco_run.log"
    shell:
        """
        busco -i {input} -o {output} \
            -l {params.busco_dir} -m genome --augustus \
            -c {threads} &> {log}
        """

rule busco_gather:
    input:
        "busco_output"
    output:
        "augustus_model/busco_genes.faa"
    params:
        lineage = config.optional.lineage,
        gene_type = "single"
    log:
        "logs/busco_gather.log"
    script:
        "scripts/busco_complete_aa.py &> {log}"

rule clustering_busco_genes:
    input:
        "augustus_model/busco_genes.faa"
    output:
        "augustus_model/cdhit.lst"
    conda:
        "envs/busco.yaml"
    log:
        "logs/clustering_busco_genes.log"
    shell:
        """
        cd-hit -o complete_buscos.cdhit -c 0.8 -i {input} -p 1 -d 0 -T 4 -M 48000 &> {log}
        grep ">" complete_buscos.cdhit | cut -f2 -d">" | cut -f1 > {output}
        """

rule concatenate_gff:
    input:
        "augustus_model/cdhit.lst",
        "busco_output"
    output:
        directory("augustus_model/busco_genes")
    conda:
        "envs/busco.yaml"
    params:
        lineage = config.optional.lineage,
        gene_type = "single"
    log:
        "logs/concatenate_gff.log"
    script:
        "scripts/concatenate_GFF.py &> {log}"

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
    log:
        "logs/gtf2genbank.log"
    shell:
        """
        gff2gbSmallDNA.pl {input.gff} {input.genome} {params.flanking_region} {output} &> {log}
        """

rule generate_subsets:
    input:
        "augustus_model/busco_genes/busco_genes.gb"
    output:
        directory("augustus_model/subsets")
    params:
        subset = calculate_subset_size("augustus_model/busco_genes/busco_genes.gb", config.optional.test_size),
        seed = 123
    log:
        "logs/generate_subsets.log"
    script:
        "scripts/generate_subsets.py &> {log}"

rule new_species:
    input:
        "augustus_model/subsets/{subset}.gb"
    output:
        touch("augustus_model/subsets/{subset}_new_species.done")
    conda:
        "envs/augustus.yaml"
    log:
        "logs/{subset}_new_species.log"
    shell:
        "new_species.pl --species={wildcards.subset} &> {log}"

rule initial_etraining:
    input:
        gb = "augustus_model/subsets/{subset}.gb",
        new_species = "augustus_model/subsets/{subset}_new_species.done"
    output:
        "augustus_model/subsets/{subset}_etrain.out"
    conda:
        "envs/augustus.yaml"
    log:
        "logs/{subset}_initial_etraining.log"
    shell:
        "etraining --species={wildcards.subset} {input.gb} &> {log}"

rule identify_bad_genes:
    input:
        "augustus_model/subsets/{subset}_etrain.out"
    output:
        "augustus_model/subsets/{subset}_bad.lst"
    log:
        "logs/{subset}_identify_bad_genes.log"
    shell:
        "grep 'in sequence' {input} | cut -f7 -d' ' | sed s/://g | sort -u > {output} &> {log}"

rule filter_genes:
    input:
        bad_list = "augustus_model/subsets/{subset}_bad.lst",
        gb = "augustus_model/subsets/{subset}.gb"
    output:
        "augustus_model/subsets/{subset}.f.gb"
    conda:
        "envs/augustus.yaml"
    log:
        "logs/{subset}_filter_genes.log"
    shell:
        "filterGenes.pl {input.bad_list} {input.gb} > {output} &> {log}"

rule retrain:
    input:
        "augustus_model/subsets/{subset}.f.gb"
    output:
        "augustus_model/subsets/{subset}_etrain.f.out"
    conda:
        "envs/augustus.yaml"
    log:
        "logs/{subset}_retrain.log"
    shell:
        "etraining --species={wildcards.subset} {input} &> {log}"

rule extract_stop_codon_freq:
    input:
        "augustus_model/subsets/{subset}_etrain.f.out"
    output:
        "augustus_model/subsets/{subset}_SC_freq.txt"
    log:
        "logs/{subset}_extract_stop_codon_freq.log"
    shell:
        "tail -6 {input} | head -3 > {output} &> {log}"

rule modify_stop_codon_freq:
    input:
        "augustus_model/subsets/{subset}_SC_freq.txt"
    output:
        touch("augustus_model/subsets/{subset}_SC_freq_modified.done")
    params:
        utilities = config.optional.utilities
    conda:
        "envs/augustus.yaml"
    log:
        "logs/{subset}_modify_stop_codon_freq.log"
    shell:
        "python3 {params.utilities}/AUGUSTUS/modify_SC_freq.py {input} {wildcards.subset} &> {log}"

rule run_augustus:
    input:
        chr19 = config.mandatory.chr19,
        sc_freq_modified = "augustus_model/subsets/{subset}_SC_freq_modified.done"
    output:
        "augustus_model/subsets/{subset}.gtf"
    conda:
        "envs/augustus.yaml"
    log:
        "logs/{subset}_run_augustus.log"
    shell:
        "augustus --species={wildcards.subset} {input.chr19} --protein=off > {output} &> {log}"