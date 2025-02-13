import glob
import attrmap as ap

configfile: os.path.join(workflow.basedir, "config.yaml")
config = ap.AttrMap(config)

include: os.path.join("rules","isoannot.smk")

include: os.path.join("rules","ab_initio.smk")

rule all:
    input:
        "augustus_model/subsets/{subset}.gtf"
