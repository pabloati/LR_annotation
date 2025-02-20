'''
Script to modify th cfg file of a species to correct the frecuency of the stop
codons
'''
import os

original_path = os.environ['PATH']

def main():
    AUGUSTUS_DIR = os.path.join(os.environ.get("AUGUSTUS_CONFIG_PATH"),"species")
    d_freq = dict()
    species = snakemake.params.name
    # Get the frequency of the Stop codons from the etrain output
    with open(snakemake.input[0], "r") as f_in:
        for linea in f_in:
            l_l = linea.split()
            codon = l_l[0].strip(":")
            freq = l_l[-1].strip("()")
            d_freq[codon] = freq
    
    # Path to the config file
    config_file = os.path.join(os.environ.get("AUGUSTUS_CONFIG_PATH"),
                               "species", species, f"{species}_parameters.cfg")
    config_file_modificado = config_file + ".tmp"
    f_out = open(config_file_modificado, "w")
    # Modify the config file
    with open(config_file, "r") as f_in:
        for linea in f_in:
            if linea.startswith("/Constant/amberprob"):
                linea = linea.replace("0.33", d_freq["tag"])
            elif linea.startswith("/Constant/ochreprob"):
                linea = linea.replace("0.33", d_freq["taa"])
            elif linea.startswith("/Constant/opalprob"):
                linea = linea.replace("0.34", d_freq["tga"])
            f_out.write(linea)
    f_out.close()
    # Replace the ild config file for the new one
    os.rename(config_file_modificado, config_file)
    with open(snakemake.output[0], "w") as f_out:
        f_out.write(f"{species} config file modified\n")
    f_out.close()

if __name__=="__main__":
    main()