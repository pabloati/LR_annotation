

rule generate_proteome:
    input:
        annot = os.path.join(dir.out.evidence_driven,"{group}_prediction_renamed.gtf"), # TODO: generate function to get the annotation depending on the method used
        reference = config.required.genome,
    output:
        os.path.join(dir.out.evidence_driven,"{group}_prediction_renamed.aa")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "generate_proteome.log")
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    shell:
        """
        getAnnoFasta.pl {input.annot} --seqfile={input.reference}
        """"

# OMARk results

rule omamer:
    input:
        proteome = os.path.join(dir.out.evidence_driven,"{group}_prediction_renamed.aa"),
        omark_db = os.path.join(dir.tools_db,f"{config.qc.omark_db}.h5")
    output:
        os.path.join(dir.out.qc_omark,"{group}","{group}_res.omamer")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "omamer_{group}.log")
    conda:
        os.path.join(dir.envs, "omark.yaml")
    shell:
        """
        omamer --db {input.omark_db} --query {input.proteome} --out {output} &> {log}
        """

rule omark:
    input:
        omamer = os.path.join(dir.out.qc_omar,"{group}","{group}.omamer"), # TODO: generate function to get the annotation depending on the method used
        omark_db = os.path.join(dir.tools_db,f"{config.qc.omark_db}.h5")
    output:
        os.path.join(dir.out.qc_omark,"{group}","{group}.pdf")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "omark_{group}.log")
    conda:
        os.path.join(dir.envs, "omark.yaml")
    shell:
        """
        omark -f {input.omamer} -d {input.omark_db} -o $(dirname {output}) &> {log}
        """

rule busco_qc:
    input:
        proteome = os.path.join(dir.out.evidence_driven,"{group}_prediction.aa"),        
    output:
        directory(os.path.join(dir.out.qc_busco,"{group}")),
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.busco.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.medium.time
    params:
        busco_dir = dir.out.busco,
        lineage = config.ab_initio.lineage,
    log: 
        os.path.join(dir.logs, "busco_qc_{group}.log")
    conda:
        os.path.join(dir.envs, "busco.yaml")
    shell:
        """
        busco -i {params.busco_dir} -o {output} -l {params.lineage} \
            -m proteins -c {params.busco_threads} --force --download_path {dir.busco_dir} &> {log}
        """