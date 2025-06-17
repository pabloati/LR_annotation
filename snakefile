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
        expand(os.path.join(dir.out.evidence_driven,"{group}_clean_prediction.gff"),group=groups),
        expand(os.path.join(dir.out.qc_omark,"{group}","{group}.pdf"),group=groups),
        expand(os.path.join(dir.out.qc_busco,"{group}"),group=groups),
        expand(os.path.join(dir.out.qc_agat,"{group}","{group}_stats.txt"),group=groups),
