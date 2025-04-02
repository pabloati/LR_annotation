
rule busco_qc:
    input:
        ""
    output:
        directory(dir.out.qc_busco)
    params:
        busco_dir = dir.out.busco,
        busco_out = dir.out.qc_busco,
        busco_db = config.busco.db,
        busco_mode = config.busco.mode,
        busco_threads = config.busco.threads
    log: 
        os.path.join(dir.logs, "busco_qc.log")
    conda:
        os.path.join(dir.envs, "busco.yaml")
    shell:
        """
        busco -i {params.busco_dir} -o {params.busco_out} -l {params.busco_db} \
            -m {params.busco_mode} -c {params.busco_threads} --force --download_path {dir.tools_busco}
        """

rule generate_proteome:
    input:
        annot = get_final_annotation(config), # TODO: generate function to get the annotation depending on the method used
        reference = config.required.genome,
    output:
        dir.out.final_results,
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
rule omark:
    input:
        get_final_annotation(config), # TODO: generate function to get the annotation depending on the method used