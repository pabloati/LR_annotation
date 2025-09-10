import glob
import attrmap as ap

config = ap.AttrMap(config)

localrules: all, install_tama, install_sqanti
# Setup rules
include: os.path.join("rules","setup","directories.smk")
include: os.path.join("rules","setup","installations.smk")
include: os.path.join("rules","setup","functions.smk")

samples,groups = get_samples_and_groups(config.required.samples_setup)
genome_name = get_genome_name(config.required.genome)


include: os.path.join("rules","transcript_modelling.smk")

include: os.path.join("rules","ab_initio.smk")

include: os.path.join("rules","evidence_driven.smk")

include: os.path.join("rules","quality_control.smk")

rule all:
    input:
        expand(os.path.join(dir.out.evidence_driven,"{sample}_clean_prediction.gff"),sample=samples),
        expand(os.path.join(dir.out.qc_omark,"{sample}","{sample}.pdf"),sample=samples),
        expand(os.path.join(dir.out.qc_busco,"{sample}"),sample=samples),
        expand(os.path.join(dir.out.qc_agat,"{sample}","{sample}_stats.txt"),sample=samples),
