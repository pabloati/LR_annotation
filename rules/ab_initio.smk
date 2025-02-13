# Snakefile for ab initio gene prediction
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
        genome = "data/genome.fasta"
    output:
        directory(dir.out.ab_initio.busco)
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
        dir.out.ab_initio.busco
    output:
        genes=os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.faa")
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
        directory(dir.out.ab_initio.augustus_model)
    conda:
        os.path.join(dir.env,"busco.yaml")
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
        gff = os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.gff")
    output:
        gen_bank = os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.gb")
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
        gen_bank = os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.gb")
    output:
        gen_bank = os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.subset.gb")
    params:
        subset = min(config.optional.test_size,
                     calculate_gene_number("augustus_model/busco_genes/busco_genes.gb")),
        seed = 123
    log:
        "logs/generate_subsets.log"
    script:
        "scripts/generate_subsets.py &> {log}"

rule new_species:
    input:
        gen_bank = os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.subset.gb")
    output:
        touch(os.path.join(dir.out.ab_initio.augustus_model,"new_species.done"))
    conda:
        "envs/augustus.yaml"
    log:
        "logs/new_species.log"
    shell:
        "new_species.pl --species=new_species &> {log}"

rule initial_etraining:
    input:
        gb = os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.subset.gb")
        new_species = os.path.join(dir.out.ab_initio.augustus_model,"new_species.done")
    output:
        training = os.path.join(dir.out.ab_initio.augustus_training,"etrain.out")
    conda:
        "envs/augustus.yaml"
    log:
        "logs/initial_etraining.log"
    shell:
        "etraining --species=new_species {input.gb} &> {log}"

rule identify_bad_genes:
    input:
        training = os.path.join(dir.out.ab_initio.augustus_training,"etrain.out")
    output:
        bad = os.path.join(dir.out.ab_initio.augustus_training,"bad.lst")
    log:
        "logs/identify_bad_genes.log"
    shell:
        "grep 'in sequence' {input} | cut -f7 -d' ' | sed s/://g | sort -u > {output} &> {log}"

rule filter_genes:
    input:
        bad = os.path.join(dir.out.ab_initio.augustus_training,"bad.lst")
        gb = os.path.join(dir.out.ab_initio.augustus_model,"busco_genes.subset.gb")
    output:
        filt = os.path.join(dir.out.ab_initio.augustus_training,"filtered.gb")
    conda:
        "envs/augustus.yaml"
    log:
        "logs/filter_genes.log"
    shell:
        "filterGenes.pl {input.bad_list} {input.gb} > {output} &> {log}"

rule retrain:
    input:
        bad = os.path.join(dir.out.ab_initio.augustus_training,"filtered.gb")
    output:
        train = os.path.join(dir.out.ab_initio.augustus_training,"etrain_filtered.out")
    conda:
        "envs/augustus.yaml"
    log:
        "logs/retrain.log"
    shell:
        "etraining --species=new_species {input} &> {log}"

rule extract_stop_codon_freq:
    input:
        train = os.path.join(dir.out.ab_initio.augustus_training,"etrain_filtered.out")
    output:
        train = os.path.join(dir.out.ab_initio.augustus_training,"SC_freq.txt")
    log:
        "logs/extract_stop_codon_freq.log"
    shell:
        "tail -6 {input} | head -3 > {output} &> {log}"

# TODO: Rewrite this rule to be adapted to snakemake nature. It creates a new file and then substitutes the frequency one.
rule modify_stop_codon_freq:
    input:
        train = os.path.join(dir.out.ab_initio.augustus_training,"SC_freq.txt")
    output:
        mod = os.path.join(dir.out.ab_initio.augustus_training,"SC_freq_mod.done")
    params:
        utilities = config.optional.utilities
    conda:
        "envs/augustus.yaml"
    log:
        "logs/modify_stop_codon_freq.log"
    script:
        "scripts/modify_SC_freq.py"

rule run_augustus:
    input:
        chr19 = config.mandatory.chr19,
        mod = os.path.join(dir.out.ab_initio.augustus_training,"SC_freq_mod.done")
    output:
        gtf = os.path.join(dir.out.ab_inition.augustus,"ab_initio_prediction.gtf")
    conda:
        "envs/augustus.yaml"
    log:
        "logs/run_augustus.log"
    shell:
        "augustus --species=new_species {input.chr19} --protein=off > {output} &> {log}"