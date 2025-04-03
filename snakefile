import glob
import attrmap as ap

config = ap.AttrMap(config)

localrules: all, install_tama, install_sqanti
# Setup rules
include: os.path.join("rules","setup","directories.smk")
include: os.path.join("rules","setup","installations.smk")
include: os.path.join("rules","setup","functions.smk")

samples,groups = get_samples_and_groups(config.required.samples_setup)


include: os.path.join("rules","isoseq.smk")

include: os.path.join("rules","ab_initio.smk")

include: os.path.join("rules","evidence_driven.smk")

include: os.path.join("rules","quality_control.smk")

rule all:
    input:
        expand(os.path.join(dir.out.evidence_driven,"{group}_clean_prediction.gtf"),group=groups),
        expand(os.path.join(dir.out.qc_omark,"{group}","{group}.pdf"),group=groups),
        expand(os.path.join(dir.out.qc_busco,"{group}"),group=groups)
