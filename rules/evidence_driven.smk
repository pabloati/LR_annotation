
rule run_sqanti:
    input:
        isoforms = os.path.join(dir.out.isoquant,"{sample}","{sample}.transcript_models.gtf"),
        ref_gff = get_sqanti_gtf(config),
        ref_genome = config.required.genome,
        sqanti = os.path.join(dir.tools_sqanti,"sqanti_installed.done")
    output:
        classification = os.path.join(dir.out.ed_sqanti,"{sample}","{sample}_classification.txt"),
        gtf = os.path.join(dir.out.ed_sqanti,"{sample}","{sample}_corrected.cds.gtf")
    threads:
        config.resources.medium.cpus,
    conda:
        f"{dir.envs}/sqanti3.yaml"
    log:
        os.path.join(dir.logs,"run_sqanti_{sample}.log")
    resources:
        slurm_extra = f"'--qos={config.resources.medium.qos}'",
        cpus_per_task = config.resources.busco.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.medium.time
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        python {dir.tools_sqanti}/sqanti3_qc.py --isoforms {input.isoforms} --refGTF {input.ref_gff} --refFasta {input.ref_genome} \
            --dir {dir.out.ed_sqanti}/{wildcards.sample} --output {wildcards.sample} -t {threads} &> {log}
        mv {dir.out.ed_sqanti}/{wildcards.sample}/{wildcards.sample}_corrected.cds.gff3 {output.gtf}
        """

rule filter_isoforms:
    input:
        classification = os.path.join(dir.out.ed_sqanti,"{sample}","{sample}_classification.txt"),
        gtf = os.path.join(dir.out.ed_sqanti,"{sample}","{sample}_corrected.cds.gtf")
    output:
        gtf = os.path.join(dir.out.ed_sqanti,"{sample}","{sample}.filtered.gtf")
    conda:
        os.path.join(dir.envs,"sqanti3.yaml")
    log:
        os.path.join(dir.logs,"filter_sqanti_{sample}.log")
    threads:
        config.resources.small.cpus,
    params:
        json_rules = config.sqanti.json_rules 
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        python {dir.tools_sqanti}/sqanti3_filter.py rules --sqanti_class {input.classification} --filter_gtf {input.gtf} \
            -j {params.json_rules} --dir {dir.out.ed_sqanti}/{wildcards.sample} \
            --output {wildcards.sample} &> {log}
        """

rule extract_hints:
    input:
        os.path.join(dir.out.ed_sqanti,"{sample}","{sample}.filtered.gtf")
    output:
        os.path.join(dir.out.ed_hints,"{sample}","{sample}.hints.gff")
    conda:
        os.path.join(dir.envs,"busco.yaml")
    params:
        utr = config.augustus.utr
        #TODO: Perhaps add techonolgy and priority options
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        tmp_dir=$(dirname {output})/tmp
        mkdir -p $tmp_dir
        if [ {params.utr} == "True" ]; then
            grep -P "\t(CDS|exon)\t" {input} | gtf2gff.pl --printIntron --out=$tmp_dir/tmp.gff
            grep -P "\t(CDS|intron|exon)\t" $tmp_dir/tmp.gff > $tmp_dir/tmp2.gff

        else
            grep -P "\t(CDS)\t" {input} | gtf2gff.pl --printIntron --out=$tmp_dir/tmp.gff
            grep -P "\t(CDS|intron)\t" $tmp_dir/tmp.gff > $tmp_dir/tmp2.gff
        fi
        # Remove gene_id and change transcript id for grp_id
        sed -i 's/gene_id[^;]*;//g' $tmp_dir/tmp2.gff
        sed -i 's/transcript_id \\"/grp=/g' $tmp_dir/tmp2.gff
        # Add the source
        cat $tmp_dir/tmp2.gff | sed "s/\\";/;pri=1;src=PB/g" > {output}
        rm -r $tmp_dir
        """

if config.augustus.mode == "split":
    include: "split_augustus.smk"

else:
    rule augustus_hints:
        input:
            genome = config.required.genome,
            mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done"),
            gff = os.path.join(dir.out.ed_hints,"{sample}","{sample}.hints.gff")
        output:
            os.path.join(dir.out.ed_augustus,"{sample}","{sample}_prediction.gff")
        conda:
            os.path.join(dir.envs,"augustus.yaml")
        params:
            name = config.augustus.species_name,
            extcfg = config.augustus.config if config.evidence_driven.config else f"{dir.envs}/extrinsic.M.RM.PB.cfg"
        log:
            os.path.join(dir.logs,"run_augustus_ed_{sample}.log")
        resources:
            slurm_extra = f"'--qos={config.resources.big.qos}'",
            cpus_per_task = config.resources.big.cpus,
            mem = config.resources.big.mem,
            runtime =  config.resources.big.time
        shell:
            """
            augustus --species={params.name} {input.genome} --hintsfile={input.gff} \
            --extrinsicCfgFile={params.extcfg} --protein=on --codingseq=on > {output} 2> {log}
            """
