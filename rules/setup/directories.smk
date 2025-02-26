"""
Ensures the correct output paths in the final direcotry
"""
import attrmap as ap

dir = ap.AttrMap()

### BUSCO database

try:
    assert (ap.utils.to_dict(config.args)["busco_db"]) is not None
    dir.busco_db = config.args.busco_db
except (KeyError, AssertionError):
    dir.busco_db = os.path.join(workflow.basedir, "..", "databases", "busco")

### Output_location
try:
    assert config.required.outdir is not None
    dir.out.base = config.required.outdir
except (KeyError, AssertionError):
    dir.out.base = "evidence_annot"

### Workflow_dirs
dir.env = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

dir.logs = os.path.join(dir.out.base, "logs")

### Output_dirs
## Ab_initio
dir.out.ab_initio = os.path.join(dir.out.base, "ab_initio")
dir.out.ab_busco = os.path.join(dir.out.ab_initio, "busco")
dir.out.ab_augustus = os.path.join(dir.out.ab_initio, "augustus")
dir.out.ab_augustus_model = os.path.join(dir.out.ab_augustus, "model")
dir.out.ab_augustus_training = os.path.join(dir.out.ab_augustus, "training")

## IsoSeq3
dir.out.isoseq = os.path.join(dir.out.base, "isoseq")
dir.out.isoseq_lima = os.path.join(dir.out.isoseq, "lima")
dir.out.isoseq_refine = os.path.join(dir.out.isoseq, "refine")
dir.out.isoseq_cluster = os.path.join(dir.out.isoseq, "cluster")
dir.out.isoseq_mapping = os.path.join(dir.out.isoseq, "mapping")