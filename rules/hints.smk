localrules: gff2bed12, tama_setup


rule gff2bed12: # local_rule
    input:
        os.path.join(dir.tools_tama,"tama_installed.done"),
        gtf = os.path.join(dir.out.isoseq_collapsed,"{sample}","{sample}.collapsed.gff")
    output:
        os.path.join(dir.out.ed_hints,"{sample}","{sample}_filtered.bed12")
    conda:
        os.path.join(dir.env,"tama.yaml")
    threads:
        config.resources.small.cpus,
    shell:
        """
        python {dir.tools_tama}/tama_go/format_converter/tama_format_gff_to_bed12_cupcake.py {input.gtf} {output}
        """

rule tama_setup: # local_rule
    input:
        files = expand(os.path.join(dir.out.ed_hints,"{sample}","{sample}_filtered.bed12"),sample=samples),
        setup = config.required.samples_setup
    output:
        os.path.join(dir.out.ed_hints,"{group}","{group}_tama_filelist.txt")
    conda:
        os.path.join(dir.env,"sqanti3.yaml")
    threads:
        config.resources.small.cpus,
    log: 
        os.path.join(dir.logs,"tama_setup_{group}.log")
    script:
        f"{dir.scripts}/create_tama_filelist.py"


rule tama_merge: 
    input:
        filelist = os.path.join(dir.out.ed_hints,"{group}","{group}_tama_filelist.txt")
    output:
        os.path.join(dir.out.ed_hints,"{group}","{group}_merge.txt")
    conda:
        os.path.join(dir.env,"tama.yaml")
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs,"tama_merge_{group}.log")
    shell:
        """
        prefix=$(echo {output} | sed 's/_merge.txt//g')
        python {dir.tools_tama}/tama_merge.py -f {input} -p $prefix &> {log}
        """

rule tama2gtf:
    input:
        bed12 = os.path.join(dir.out.ed_hints,"{group}","{group}_merge.txt")
    output:
        os.path.join(dir.out.ed_hints,"{group}","{group}_merged.gtf")
    conda:
        os.path.join(dir.env,"tama.yaml")
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs,"tama2gtf_{group}.log")
    shell:
        """
        python {dir.tools_tama}/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py {input} {output} $>{log}
        """

rule run_sqanti:
    input:
        isoforms = os.path.join(dir.out.ed_hints,"{group}","{group}_merged.gtf"),
        ref_gff = get_sqanti_gtf(config),
        ref_genome = config.required.genome,
        sqanti = os.path.join(dir.tools_sqanti,"sqanti_installed.done")
    output:
        os.path.join(dir.out.ed_sqanti,"{group}","{group}_classification.txt"),
        os.path.join(dir.out.ed_sqanti,"{group}","{group}_corrected.gtf.cds.gff")
    threads:
        config.resources.medium.cpus,
    conda:
        f"{dir.env}/sqanti3.yaml"
    log:
        os.path.join(dir.logs,"run_sqanti_{group}.log")
    resources:
        slurm_extra = f"'--qos={config.resources.medium.qos}'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.medium.mem,
        runtime =  config.resources.medium.time
    shell:
        """
        python {dir.tools_sqanti}/sqanti3_qc.py {input.isoforms} {input.ref_gff} {input.ref_genome} \
            --dir {dir.out.ed_sqanti}/{wildcards.group} --output {wildcards.group} -t {threads} &> {log}
        """

rule filter_isoforms:
    input:
        classification = os.path.join(dir.out.ed_sqanti,"{group}","{group}_classification.txt"),
        gtf = os.path.join(dir.out.ed_sqanti,"{group}","{group}_corrected.gtf.cds.gff")
    output:
        classification = os.path.join(dir.out.ed_sqanti,"{group}","{group}_classification.filt.txt"),
        gtf = os.path.join(dir.out.ed_sqanti,"{group}","{group}_corrected.cds.filt.gtf")
    conda:
        os.path.join(dir.env,"sqanti3.yaml")
    log:
        os.path.join(dir.logs,"filter_sqanti_{group}.log")
    threads:
        config.resources.small.cpus,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        Rscript {dir.scripts}/filterClassification.R {input.classification} {input.gtf} {output.classification} {output.gtf} &> {log}
        """

rule extract_hints:
    input:
        os.path.join(dir.out.ed_sqanti,"{group}","{group}_corrected.cds.filt.gtf")
    output:
        os.path.join(dir.out.ed_hints,"{group}","{group}.hints.gff")
    conda:
        os.path.join(dir.env,"busco.yaml")
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
            gff = os.path.join(dir.out.ed_hints,"{group}","{group}.hints.gff")
        output:
            gtf = os.path.join(dir.out.evidence_driven,"{group}_prediction.gtf")
        conda:
            os.path.join(dir.env,"augustus.yaml")
        params:
            name = config.optional.species_name,
            extcfg = config.augustus.config if config.evidence_driven.config else f"{dir.envs}/extrinsic.M.RM.PB.cfg"
        log:
            os.path.join(dir.logs,"run_augustus_{group}.log")
        resources:
            slurm_extra = f"'--qos={config.resources.big.qos}'",
            cpus_per_task = config.resources.big.cpus,
            mem = config.resources.big.mem,
            runtime =  config.resources.big.time
        shell:
            """
            augustus --species={params.name} {input.genome} --hintsfile={input.gff} \
            --extrinsicCfgFile={params.extcfg} --protein=on --codingseq=on > {output}
            """
