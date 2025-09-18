import os

def get_sqanti_gtf(config):
    if config.augustus.prediction == "ab_initio":
        return os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
    elif config.augustus.prediction == "evidence_driven":
        return config.augustus.reference_gtf

def get_sample_name(file):
    sample = os.path.splitext(os.path.basename(file))[0]
    filetype = os.path.splitext(file)[1]
    return sample, filetype

def get_pbmm2_input(filetype,config):
    if filetype == "fastq":
        return config.required.input
    else:
        return os.path.join(dir.out.isoseq,f"{sample}.fastq")

def get_chromosomes(file):
    import subprocess
    result = subprocess.run(["grep", ">", file], capture_output=True, text=True)
    output_list = [line.split(' ')[0].replace('>', '') for line in result.stdout.strip().split('\n')]
    return output_list

def get_genome_name(file):
    return  os.path.splitext(os.path.basename(file))[0]

# def get_final_annotation(config):
#     if config.augustus.prediction == "ab_initio":
#         return os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
#     elif config.augustus.prediction == "evidence_driven":
#         return os.path.join(dir.out.evidence_driven,"{group}_prediction.gtf").
