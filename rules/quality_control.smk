
rule generate_proteome:
    input:
        annot = os.path.join(dir.out.ed_augustus,"{group}_prediction_renamed.gtf"), 
        reference = config.required.genome,
    output:
        os.path.join(dir.out.ed_augustus,"{group}_prediction_renamed.aa")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "generate_proteome_{group}.log")
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    shell:
        """
        getAnnoFasta.pl {input.annot} --seqfile={input.reference}
        """

# OMARk results

rule omamer:
    input:
        proteome = os.path.join(dir.out.ed_augustus,"{group}_prediction_renamed.aa"),
        omark_db = os.path.join(dir.tools_omark,f"{config.qc.omark_db}.h5")
    output:
        os.path.join(dir.out.qc_omark,"{group}","{group}.omamer")
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
        omamer = os.path.join(dir.out.qc_omark,"{group}","{group}.omamer"),
        omark_db = os.path.join(dir.tools_omark,f"{config.qc.omark_db}.h5")
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
    threads:
        config.resources.small.cpus,
    shell:
        """
        omark -f {input.omamer} -d {input.omark_db} -o $(dirname {output}) &> {log}
        """

rule busco_qc:
    input:
        proteome = os.path.join(dir.out.ed_augustus,"{group}_prediction_renamed.aa"),        
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
    threads:
        config.resources.busco.cpus,
    shell:
        """
        busco -i {params.busco_dir} -o {output} -l {params.lineage} \
            -m proteins -c {threads} --force --download_path {dir.busco_dir} &> {log}
        """