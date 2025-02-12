import glob

configfile: os.path.join(workflow.basedir, "config.yaml")

include: os.path.join("rules","isoannot.smk")

include: os.path.join("rules","ab_initio.smk")

rule all:
    input:
        