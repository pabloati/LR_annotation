# SQANTI-evidence

**Long-Read Evidence-Driven Structural Annotation Pipeline**  
A Snakemake workflow to produce structural genome annotations leveraging long-read sequencing data.

---

## Table of Contents

1. [Overview](#overview)  
2. [Features](#features)  
3. [Requirements](#requirements)  
4. [Installation](#installation)  
5. [Usage](#usage)  
6. [Configuration](#configuration)  
7. [Pipeline Workflow](#pipeline-workflow)  
8. [Scripts & Rules](#scripts--rules)  
9. [Output](#output)  
10. [Examples](#examples)  
11. [License & Citations](#license--citations)  
12. [Contact / Support](#contact--support)

---

## Overview

This repository implements a **Snakemake** pipeline (with auxiliary scripts) to generate structural genome annotations guided by long-read sequencing data (e.g. PacBio, Oxford Nanopore).  
It aims to produce high-quality annotations by combining transcript evidence from long reads with conventional annotation strategies. The main structure of the pipeline and use of the long-read transcriptomics is derived from [this paper](https://genome.cshlp.org/content/early/2025/03/04/gr.279864.124).

---

## Features

- Modular pipeline built with Snakemake  
- Integration of long-read data to inform exon/intron boundaries  
- Flexible configuration for different organisms & datasets  
- Support for cluster execution (e.g. SLURM)  
- Scripts to assist in annotation processing and QC  

---

## Requirements

- **Snakemake** (version  >= X.X)  
- **Python** (>= 3.8) + dependencies  
- Linux / Unix environment  
- Long-read RNA (or cDNA) sequencing aligned data (BAM or SAM)  
- Reference genome (FASTA)  
- (Optional) Annotation hints / protein / transcript evidence  

You’ll find an `envs/` folder for environment / dependency configurations.

---

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/pabloati/LR_annotation.git
   cd LR_annotation
````

2. Create and activate a conda / mamba environment (if using):

   ```bash
   conda env create -f envs/env.yaml
   conda activate <env_name>
   ```

   (You may have multiple environments defined under `envs/`, inspect and choose the appropriate one.)

3. Install any extra Python packages not handled by the environment file:

   ```bash
   pip install -r requirements.txt
   ```

   *(If `requirements.txt` does not exist, you can generate one from the environment.)*

---

## Usage

From the root directory, run:

```bash
snakemake --profile profile_slurm.yaml  # or another execution mode
```

Or for local (no cluster) execution:

```bash
snakemake -j <num_threads>
```

You can also target specific rules or outputs:

```bash
snakemake path/to/output.gff3
```

---

## Configuration

The behaviour of the pipeline is controlled via:

* `config.yaml` — main configuration file (genome paths, sample IDs, parameters)
* `profile_slurm.yaml` — parameters and settings for SLURM (if using cluster)

Edit `config.yaml` to point to your reference genome, aligned reads, and other evidence files.

---

## Pipeline Workflow

Rough outline of the major steps / rules (in `rules/`):

1. Preprocessing of reads / alignments
2. Transcript feature extraction
3. Long-read informed exon/intron boundary refinement
4. Evidence merging with other annotation sources
5. Final structural annotation (e.g. GFF3 output)
6. QC and filtering steps

Refer to the individual rule files in `rules/` for detailed logic.

---

## Scripts & Rules

* `scripts/` — utility scripts used by the workflow (e.g. parsing, filtering)
* `snakefile` — main workflow entry
* `rules/` — subrules modularizing steps
* `lr_annot.py` — core Python module / driver (if used in pipeline)

You can read through them to see custom parameters, function calls, and expected behavior.

---

## Output

Typical outputs include:

* GFF3 / GTF annotated structural models
* Transcript / exon / intron files
* QC reports
* Intermediate alignment / feature files

Output paths and filenames are configurable via `config.yaml`.

---

## Examples

*(You may want to include a small example or test dataset to demonstrate pipeline execution. If you have one, mention it here. E.g.):*

* `example/` — folder with toy genome + reads, config, and expected outputs

* Usage:

  ```bash
  cd example
  snakemake -j 4
  ```

* Compare output GFF3 with expected reference.

If you don’t have an example, you could add one in future to help users.

---

## License & Citations

State your license (e.g. MIT, GPL, etc.) here.
Also include citations to relevant tools or papers used in this pipeline.

```text
MIT License
(c) 2025 Pablo A. Oti (or your name)
```

Please cite this repository as:

> Oti, P. (2025). **LR_annotation**: Long-Read Guided Structural Annotation Pipeline. GitHub. [https://github.com/pabloati/LR_annotation](https://github.com/pabloati/LR_annotation)

---

## Contact / Support

For questions or issues, open an **Issue** on GitHub.
You can also reach me at: *[your_email@example.com](mailto:your_email@example.com)*.

---

## Future Improvements

You might want to add:

* A Docker or Singularity container for reproducibility
* Automated tests / CI
* Support for additional evidence types
* Visualization modules
* More extensive examples & documentation

---

That’s the baseline README. If you send me specific details (license, example dataset, more description of what rules do), I can further customize it.

