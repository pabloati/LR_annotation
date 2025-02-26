import glob
import attrmap as ap

configfile: os.path.join(workflow.basedir, "config.yaml")
config = ap.AttrMap(config)

# Setup rules
include: os.path.join("rules","setup","directories.smk")


include: os.path.join("rules","isoseq.smk")

include: os.path.join("rules","ab_initio.smk")

rule all:
    input:
        os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf"),
        expand(os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.fastq"),sample=samples)
