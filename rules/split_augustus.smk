chromosomes=get_chromosomes(config.required.genome)
rule split_fasta:
    input:
        fasta = config.required.genome
    output:
        touch(os.path.join(dir.tools_reference,"fasta_split.done"))
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    conda:
        f"{dir.env}/sqanti3.yaml"
    threads:
        config.resources.small.cpus
    log:
        os.path.join(dir.logs,"split_fasta.log")
    script:
        os.path.join(dir.scripts,"splitfasta.py")

rule ed_augusuts_per_chromosome:
    input:
        os.path.join(dir.tools_reference,"fasta_split.done"),
        mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done"),
        gff = os.path.join(dir.out.ed_hints,"{group}","{group}.hints.gff")
    output:
        os.path.join(dir.out.evidence_driven,"augustus","{group}","{chromosome}.prediction.gff")
    conda:
        f"{dir.env}/augustus.yaml"
    params:
        name = config.optional.species_name,
        extcfg = f"{dir.envs}/extrinsic.M.RM.PB.cfg"
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    threads:
        config.resources.small.cpus
    log:
        os.path.join(dir.logs,"run_augustus_{group}_{chromosome}.log")
    shell:
        """
        chromosome={dir.tools_reference}/{wildcards.chromosome}.fasta
        augustus --species={params.name} $chromosome --hintsfile={input.gff} \
        --extrinsicCfgFile={params.extcfg} --protein=on --codingseq=on > {output} 
        """
rule merge_ed_predicitons:
    input:
        expand(os.path.join(dir.out.evidence_driven,"augustus","{{group}}","{chromosome}.prediction.gff"),chromosome=chromosomes)
    output:
        gtf = os.path.join(dir.out.evidence_driven,"{group}_prediction.gtf")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    threads:
        config.resources.small.cpus
    shell:
        """
        for file in {input} ; do
            if [[ $((grep -v -c "#" $file)) -gt 0 ]]; then
                cat $file >> {output}
            fi
        done
        """

rule ab_augustus_per_chromosome:
    input:
        os.path.join(dir.tools_reference,"fasta_split.done"),
        mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done")
    output:
        os.path.join(dir.out.ab_augustus,"split","{chromosome}.prediction.gff")
    conda:
        f"{dir.env}/augustus.yaml"
    params:
        name = config.optional.species_name,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    threads:
        config.resources.small.cpus
    log:
        os.path.join(dir.logs,"run_augustus_{chromosome}.log")
    shell:
        """
        chromosome={dir.tools_reference}/{wildcards.chromosome}.fasta
        augustus --species={params.name} $chromosome --protein=on --codingseq=on > {output} 
        """

rule merge_ab_predictions:
    input:
        expand(os.path.join(dir.out.ab_augustus,"split","{chromosome}.prediction.gff"),chromosome=chromosomes)
    output:
        gtf = os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    threads:
        config.resources.small.cpus
    shell:
        """
        cat {input} | grep -v "#" > {output}
        """

