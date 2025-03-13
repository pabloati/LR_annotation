import os

rule run_sqanti:
    input:
        isoforms = os.path.join(dir.out.isoseq_collapsed,"{sample}","{sample}.collapsed.gff"),
        ab_initio_gff = get_sqanti_gtf(config),
        ref_genome = config.required.genome,
    output:
        os.path.join(dir.out.sqanti,"{sample}","{sample}_classification.txt"),
        os.path.join(dir.out.sqanti,"{sample}","{sample}_corrected.gtf")
    conda:
        f"{dir.env}/sqanti3.yaml"
    resources:
        slurm_extra = f"'--qos={config.resources.medium.qos}'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.medium.mem,
        runtime =  config.resources.medium.time
    shell:
        """
        sqanti_qc.py {input.isoforms} {input.ab_initio_gff} {input.ref_genome} \
            --dir $(basename {output}) --prefix {wildcards.sample} \
            --skipORF
        """

rule filter_isoforms:
    input:
        classification = os.path.join(dir.out.sqanti,"{sample}","{sample}_classification.txt"),
        gtf = os.path.join(dir.out.sqanti,"{sample}","{sample}_corrected.gtf")
    output:
        classification = os.path.join(dir.out.sqanti,"{sample}_classification.filt.txt"),
        gtf = os.path.join(dir.out.evidence_driven,"{sample}","{sample}_filtered.gtf")
    conda:
        os.path.join(dir.env,"sqanti3.yaml")
    log:
        os.path.join(dir.logs,"filter_sqanti_{sample}.log")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    script:
        os.path.join(dir.scripts,"filterClassifications.R")

rule gff2bed12:
    input:
        os.path.join(dir.tools_tama,"tama_installed.done"),
        gtf = os.path.join(dir.out.evidence_driven,"{sample}","{sample}_filtered.gtf")
    output:
        os.path.join(dir.out.evidence_driven,"{sample}","{sample}_filtered.bed12")
    conda:
        os.path.join(dir.env,"tama.yaml")
    shell:
        """
        python {dir.tools_tama}/tama_go/format_converter/tama_format_gff_to_bed12_cupcake.py {input.gtf} {output}
        """

rule tama_setup:
    input:
        files = expand(os.path.join(dir.out.evidence_driven,"{sample}","{sample}_filtered.bed12"),sample=samples),
        setup = config.required.samples_setup
    output:
        os.path.join(dir.out.evidence_driven,"{group}_tama_filelist.txt")
    conda:
        os.path.join(dir.env,"tama.yaml")
    script:
        os.path.join(dir.scripts,"tama_setup.py")

rule tama_merge:
    input:
        filelist = os.path.join(dir.out.evidence_driven,"{group}_tama_filelist.txt")
    output:
        os.path.join(dir.out.evidence_driven,"{group}_merged.bed12")
    shell:
        """
        python {dir.tools_tama}/tama_go/format_converter/tama_merge_bed12.py {input} {output}
        """

rule tama2gtf:
    input:
        bed12 = os.path.join(dir.out.evidence_driven,"{group}_merged.bed12")
    output:
        os.path.join(dir.out.evidence_driven,"{group}_merged.gtf")
    conda:
        os.path.join(dir.env,"tama.yaml")
    shell:
        """
        python {dir.tools_tama}/tama_go/format_converter/tama_format_bed12_to_gtf.py {input} {output}
        """

rule ORF_prediction:
    input:
        gtf = os.path.join(dir.out.evidence_driven,"{group}_merged.gtf")
    output:
        os.path.join(dir.out.evidence_driven,"{group}_merged_cds.gtf")
    conda:
        os.path.join(dir.env,"tama.yaml")
    shell:
        """
        GeneMarkS -m {input} -o {output}
        """

rule extract_hints:
    input:
        os.path.join(dir.out.evidence_driven,"{group}","{group}_merged_cds.gtf")
    output:
        os.path.join(dir.out.evidence_driven,"{group}","{group}.hints.gff")
    conda:
        os.path.join(dir.env,"busco.yaml")
    params:
        utr = config.hints.utr
        #TODO: Perhaps add techonolgy and priority options
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        if [ {params.utr} == "True" ]; then
            grep -P "\t(CDS|exon)\t" {input} | gtf2gff.pl --printIntron --out=tmp.gff
            grep -P "\t(CDS|intron|exon)\t" tmp.gff > tmp2.gff

        else
            grep -P "\t(CDS)\t" $gff | gtf2gff.pl --printIntron --out=tmp.gff
            grep -P "\t(CDS|intron)\t" tmp.gff > tmp2.gff
        fi
        # Remove gene_id and change transcript id for grp_id
        sed -i 's/gene_id[^;]*;//g' tmp2.gff
        sed -i 's/transcript_id \"/grp=/g' tmp2.gff

        # change the trancript id for the source
        cat tmp2.gff | sed "s/\";$/;pri=1;src=PB/g" > {output}
        rm tmp.gff tmp2.gff
        """

rule augustus_hints:
    input:
        genome = config.required.genome,
        mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done"),
        gff = os.path.join(dir.out.evidence_driven,"{group}","{group}.hints.gff")
    output:
        gtf = os.path.join(dir.out.evidence_driven,"{group}","evidence_driven_prediction.gtf")
    conda:
        os.path.join(dir.env,"augustus.yaml")
    params:
        name = config.optional.species_name
    log:
        os.path.join(dir.logs,"run_augustus_{group}.log")
    resources:
        slurm_extra = f"'--qos={config.resources.big.qos}'",
        cpus_per_task = config.resources.big.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    shell:
        "augustus --species={params.name} {input.genome} --hintsfile={input.gff} --protein=off > {output} &> {log}"
