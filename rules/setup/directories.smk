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
    assert (ap.utils.to_dict(config.args)["required"]["output"]) is not None
    dir.out.base = config.args.output
except (KeyError, AssertionError):
    dir.out.base = "evidence_annot"

### Workflow_dirs
dir.env = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

### Output_dirs
## Ab_initio
dir.out.ab_initio = os.path.join(dir.out.base, "ab_initio")
dir.out.ab_initio.busco = os.path.join(dir.out.ab_initio, "busco")
dir.out.ab_initio.augustus = os.path.join(dir.out.ab_initio, "augustus")
dir.out.ab_initio.augustus_model = os.path.join(dir.out.ab_initio.augustus, "model")
dir.out.ab_initio.augustus_training = os.path.join(dir.out.ab_initio.augustus, "training")
## IsoSeq3
