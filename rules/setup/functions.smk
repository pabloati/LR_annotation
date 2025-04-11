import os

def get_sqanti_gtf(config):
    if config.augustus.prediction == "ab_initio":
        return os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
    elif config.augustus.prediction == "evidence_driven":
        return config.augustus.reference_gtf

def get_samples_and_groups(file):
    samples = []
    groups = []
    with open(file) as f:
        i=0
        for line in f:
            if i == 0:
                i+=1
                continue
            else:
                samples.append(line.split("\t")[0])
                if line.split("\t")[1] not in groups:
                    groups.append(line.split("\t")[1])
    return samples, groups

def get_chromosomes(file):
    import subprocess
    result = subprocess.run(["grep", ">", file], capture_output=True, text=True)
    output_list = [line.replace('>', '') for line in result.stdout.strip().split('\n')]
    return output_list

def get_genome_name(file):
    return  os.path.splitext(os.path.basename(file))[0]

# def get_final_annotation(config):
#     if config.augustus.prediction == "ab_initio":
#         return os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
#     elif config.augustus.prediction == "evidence_driven":
#         return os.path.join(dir.out.evidence_driven,"{group}_prediction.gtf").

def get_tama_input(approach):
    if approach == "isoseq3":
        return os.path.join(dir.out.isoseq_collapsed,"{sample}","{sample}.collapsed.gff")
    elif approach == "isoquant":
        return os.path.join(dir.out.isoquant,"{sample}","{sample}.transcript_models.gtf")

def get_input_type(filename):
    if filename.endswith(".bam"):
        return "bam"
    elif filename.endswith(".fastq"):
        return "fasta"
    elif filename.endswith(".fasta") or filename.endswith(".fa"):
        return "fasta"
    else:
        raise ValueError("Unsupported file type. Please provide a .bam, .fastq, or .fasta file.")
    
def get_sqanti_input(approach,size):
    if size == 1:
        if approach == "isoseq3":
            return os.path.join(dir.out.isoseq_collapsed,"{group}","{group}.collapsed.gff")
    elif approach == "isoquant":
            return os.path.join(dir.out.isoquant,"{group}","{group}.transcript_models.gtf")
    elif size == 2:
        return os.path.join(dir.out.ed_hints,"{group}","{group}_merged.gtf") 