
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
        mem = config.resources.medium.cpus,
        runtime =  config.resources.medium.time
    shell:
        """
        sqanti_qc.py {input.isoforms} {input.ab_initio_gff} {input.ref_genome} --dir $(basename {output}) --prefix {wildcards.sample}
        """