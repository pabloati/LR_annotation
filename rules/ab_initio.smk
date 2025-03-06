# Snakefile for ab initio gene prediction

# Setup local rules (do not require much resources)
localrules: new_species, identify_bad_genes, extract_stop_codon_freq
#TODO: move to a module
def calculate_gene_number(file_path):
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


rule busco_run:
    input:
        genome = config.required.genome
    output:
        directory(dir.out.ab_busco)
    conda:
        f"{dir.env}/busco.yaml"  
    params:
        busco_dir = config.optional.busco_downloads,
        lineage = config.optional.lineage
    resources:
        qos = config.resources.busco.qos,
        cpus_per_task = config.resources.busco.cpus,
        mem = config.resources.busco.cpus,
        runtime =  config.resources.busco.time
    threads:
        config.resources.big.cpus,
    log:
        os.path.join(dir.logs,"busco_run.log")
    shell:
        """
        busco -i {input} -o {output} \
            -l {params.lineage} -m genome --augustus \
            -c {threads} --download_path {params.busco_dir} &> {log}
        """

rule busco_gather:
    input:
        dir.out.ab_busco
    output:
        genes=os.path.join(dir.out.ab_augustus_model,"busco_genes.faa")
    params:
        lineage = config.optional.lineage,
        gene_type = "single"
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    conda:
        os.path.join(dir.env,"busco.yaml")
    log:
        os.path.join(dir.logs,"busco_gather.log")
    script:
        os.path.join(dir.scripts,"busco_complete_aa.py")

rule clustering_busco_genes:
    input:
        os.path.join(dir.out.ab_augustus_model,"busco_genes.faa")
    output:
        os.path.join(dir.out.ab_augustus_model,"cdhit.lst")
    conda:
        os.path.join(dir.env,"busco.yaml")
    log:
        os.path.join(dir.logs,"clustering_busco_genes.log")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    shell:
        """
        dir=$(dirname {output})
        cd-hit -o $dir/complete_buscos.cdhit -c 0.8 -i {input} -p 1 -d 0 -T 4 -M 48000 &> {log}
        grep ">" $dir/complete_buscos.cdhit | cut -f2 -d">" | cut -f1 > {output}
        """

rule concatenate_gff:
    input:
        gene_list = os.path.join(dir.out.ab_augustus_model,"cdhit.lst"),
        busco_path = dir.out.ab_busco
    output:
        os.path.join(dir.out.ab_augustus_model,"busco_genes.gff")
    conda:
        os.path.join(dir.env,"busco.yaml")
    params:
        lineage = config.optional.lineage,
        gene_type = "single"
    log:
        os.path.join(dir.logs,"concatenate_gff.log")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    script:
        os.path.join(dir.scripts,"concatenate_GFF.py")

rule gtf2genbank:
    input:
        genome = config.required.genome,
        gff = os.path.join(dir.out.ab_augustus_model,"busco_genes.gff")
    output:
        gen_bank = os.path.join(dir.out.ab_augustus_model,"busco_genes.gb")
    conda:
        os.path.join(dir.env,"busco.yaml")
    params:
        flanking_region = config.optional.flanking_region
    log:
        os.path.join(dir.logs,"gtf2genbank.log")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    shell:
        """
        gff2gbSmallDNA.pl {input.gff} {input.genome} {params.flanking_region} {output} &> {log}
        """

rule generate_subsets:
    input:
        gen_bank_in = os.path.join(dir.out.ab_augustus_model,"busco_genes.gb")
    output:
        gen_bank_out = os.path.join(dir.out.ab_augustus_model,"busco_genes.subset.gb")
    params:
        size = config.optional.test_size,
        seed = 123
    log:
        os.path.join(dir.logs,"generate_subset.log")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    script:
        os.path.join(dir.scripts,"generate_subset.py")

# TODO: Skip this rule if the directory of the new species exist
# TODO: Change the done path to be to a specific directory for the "check files"
rule new_species:
    input:
        gen_bank = os.path.join(dir.out.ab_augustus_model,"busco_genes.subset.gb")
    output:
        touch(os.path.join(dir.out.ab_augustus_model,f"{config.optional.species_name}.done"))
    conda:
        os.path.join(dir.env,"augustus.yaml")
    params:
        name = config.optional.species_name,
        augustus_dir = os.environ.get("AUGUSTUS_CONFIG_PATH")
    log:
        os.path.join(dir.logs,"new_species.log")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    shell:
        """
        rm -rf {params.augustus_dir}/species/{params.name}
        new_species.pl --species={params.name} &> {log}
        """

rule initial_etraining:
    input:
        gb = os.path.join(dir.out.ab_augustus_model,"busco_genes.subset.gb"),
        new_species = os.path.join(dir.out.ab_augustus_model,f"{config.optional.species_name}.done")
    output:
        training = os.path.join(dir.out.ab_augustus_training,"etrain.out")
    conda:
        os.path.join(dir.env,"augustus.yaml")
    params:
        name = config.optional.species_name
    log:
        os.path.join(dir.logs,"initial_etraining.log")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    shell:
        "etraining --species={params.name} {input.gb} &> {output}"

rule identify_bad_genes:
    input:
        training = os.path.join(dir.out.ab_augustus_training,"etrain.out")
    output:
        bad = os.path.join(dir.out.ab_augustus_training,"bad.lst")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    shell:
        "grep 'in sequence' {input} | cut -f7 -d' ' | sed s/://g | sort -u > {output}"

rule filter_genes:
    input:
        bad_list = os.path.join(dir.out.ab_augustus_training,"bad.lst"),
        gb = os.path.join(dir.out.ab_augustus_model,"busco_genes.subset.gb")
    output:
        filt = os.path.join(dir.out.ab_augustus_training,"filtered.gb")
    conda:
        os.path.join(dir.env,"augustus.yaml")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    shell:
        "filterGenes.pl {input.bad_list} {input.gb} > {output}"

rule retrain:
    input:
        bad = os.path.join(dir.out.ab_augustus_training,"filtered.gb")
    output:
        train = os.path.join(dir.out.ab_augustus_training,"etrain_filtered.out")
    conda:
        os.path.join(dir.env,"augustus.yaml")
    params:
        name = config.optional.species_name
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    shell:
        "etraining --species={params.name} {input} > {output}"

rule extract_stop_codon_freq:
    input:
        train = os.path.join(dir.out.ab_augustus_training,"etrain_filtered.out")
    output:
        train = os.path.join(dir.out.ab_augustus_training,"SC_freq.txt")
    shell:
        "tail -6 {input} | head -3 > {output}"

# TODO: Rewrite this rule to be adapted to snakemake nature. It creates a new file and then substitutes the frequency one.
# TODO: Find a way to give the shell variable AUGUSUTUS_CONFIG_PATH directly to the file
rule modify_stop_codon_freq:
    input:
        train = os.path.join(dir.out.ab_augustus_training,"SC_freq.txt")
    output:
        mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done")
    params:
        name = config.optional.species_name
    conda:
        os.path.join(dir.env,"augustus.yaml")
    log:
        os.path.join(dir.logs,"modify_stop_codon_freq.log")
    resources:
        qos = config.resources.small.qos,
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.cpus,
        runtime =  config.resources.small.time
    script:
        os.path.join(dir.scripts,"modify_SC_freq.py")

# TODO: Is there any way to increase augustus usage to >1 core?
rule run_augustus:
    input:
        genome = config.required.genome,
        mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done")
    output:
        gtf = os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
    conda:
        os.path.join(dir.env,"augustus.yaml")
    params:
        name = config.optional.species_name
    log:
        os.path.join(dir.logs,"run_augustus.log")
    resources:
        qos = config.resources.big.qos,
        cpus_per_task = config.resources.big.cpus,
        mem = config.resources.big.cpus,
        runtime =  config.resources.big.time
    shell:
        "augustus --species={params.name} {input.genome} --protein=off > {output} &> {log}"