import glob
import attrmap as ap

config = ap.AttrMap(config)

localrules: all
# Setup rules
include: os.path.join("rules","setup","directories.smk")
include: os.path.join("rules","setup","installations.smk")
include: os.path.join("rules","setup","functions.smk")

samples,groups = get_samples_and_groups(config.required.samples_setup)


include: os.path.join("rules","isoseq.smk")

include: os.path.join("rules","ab_initio.smk")

include: os.path.join("rules","hints.smk")

print(samples)
print(groups)
rule all:
    input:
        #expand(os.path.join(dir.out.evidence_driven,"{group}","evidence_driven_prediction.gtf"),group=groups)
        expand( os.path.join(dir.out.evidence_driven,"{group}","{group}.hints.gff"),group=groups)