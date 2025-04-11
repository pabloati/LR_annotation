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
