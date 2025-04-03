import os,sys

def main():
    dir=os.path.dirname(snakemake.input[0])
    samples_file = snakemake.params[0]
    brk_2_samples = {}

    with open(samples_file) as f:
        samples = f.readlines()
        for line in samples:
            line = line.strip("\n")
            brk_2_samples[line.split(",")[0]] = line.split(",")[1]
    
    for brk in os.listdir(dir):
        directory = os.path.join(dir, brk)
        if os.path.isdir(directory):
            try:
                sample_name = brk_2_samples[brk]
                path=os.path.join(dir, sample_name)
                try:
                    os.mkdir(path)
                except:
                    pass
                for file in os.listdir(directory):
                    new_name = file.replace(f"fl.{brk}", f"{sample_name}.fl")
                    os.rename(os.path.join(directory, file), 
                            os.path.join(path,new_name))
                os.rmdir(directory)
            except KeyError:
                if brk not in list(brk_2_samples.values()):
                    print(f"Sample {brk} not found in samples file")
                    sys.exit(1)
            
                


if __name__=="__main__":
    main()
    with open(snakemake.output[0], "w") as f:
        f.write("Lima has been correctly renamed\n")
    f.close()