
# rule merge_isoforms:
#     input:
#         expand(os.path.join(dir.out.isoseq_collapsed,"{sample}","{sample}.collapsed.gff"),sample=samples)
#     output:
#         os.path.join(dir.out.isoseq,"merged.gff")
#     conda:
#         f"{dir.env}/agat.yaml"
#     shell:
#         """
#         cat {input} > {output}
#         """
    
rule run_sqanti:
    input:
        isoforms = os.path.join(dir.out.isoseq_collapsed,"{sample}","{sample}.collapsed.gff"),
        ab_initio_gff = os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf"),
        ref_genome = config.required.genome,
    output:
        os.path.join(dir.out.sqanti,"{sample}","{sample}_classification.txt")
    conda:
        f"{dir.env}/sqanti3.yaml"
    resources:
        slurm_extra = f"'--qos={config.resources.medium.qos}'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.medium.mem,
        runtime =  config.resources.medium.time
    shell:
        """
        sqanti_qc.py {input.isoforms} {input.ab_initio_gff} {input.ref_genome} --dir $(basename {output}) --prefix {wildcards.sample}
        """

rule filter_isoforms:
    input:
        classification = os.path.join(dir.out.sqanti,"{sample}","{sample}_classification.txt"),
        gtf = os.path.join(dir.out.sqanti,"{sample}","{sample}_corrected.gtf.cds.gff")
    output:
        classification = os.path.join(dir.out.evidence_driven,"{sample}_classification.filt.txt"),
        gtf = os.path.join(dir.out.evidence_driven,"{sample}","{sample}filtered_cds.gtf")
    conda:
        os.path.join(dir.env,"sqanti3.yaml")
    log:
        os.path.join(dir.logs,"filter_sqanti.log")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    script:
        os.path.join(dir.scripts,"filterClassifications.R")

rule tama_setup:


rule extract_hints:
    input:
        os.path.join(dir.out.evidence_driven,"{sample}","{sample}filtered_cds.gtf")
    output:
        os.path.join(dir.out.evidence_driven,"{sample}","{sample}.hints.gff")
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
        gff = os.path.join(dir.out.hints,"{sample}","{sample}.hints.gff")
    output:
        gtf = os.path.join(dir.out.evidence_driven," {sample}","evidence_driven_prediction.gtf")
    conda:
        os.path.join(dir.env,"augustus.yaml")
    params:
        name = config.optional.species_name
    log:
        os.path.join(dir.logs,"run_augustus.log")
    resources:
        slurm_extra = f"'--qos={config.resources.big.qos}'",
        cpus_per_task = config.resources.big.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    shell:
        "augustus --species={params.name} {input.genome} --hintsfile={input.gff} --protein=off > {output} &> {log}"
