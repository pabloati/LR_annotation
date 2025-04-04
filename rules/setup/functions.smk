
def get_sqanti_gtf(config):
    if config.augustus.prediction == "ab_initio":
        return os.path.join(dir.out.ab_augustus_model,"augustus.gtf")
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