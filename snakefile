import glob
import attrmap as ap

configfile: os.path.join(workflow.basedir, "config.yaml")
config = ap.AttrMap(config)

localrules: all
# Setup rules
include: os.path.join("rules","setup","directories.smk")


include: os.path.join("rules","isoseq.smk")

include: os.path.join("rules","ab_initio.smk")

include: os.path.join("rules","hints.smk")

rule all:
    input:
        expand(os.path.join(dir.out.sqanti,"{sample}","{sample}_classification.txt"),sample=samples)
