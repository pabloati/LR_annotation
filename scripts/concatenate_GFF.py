'''
Author: Alejandro Paniagua
Modified by: Pablo Atienza
Given a file with the BUSCO complete genes that what to be incoporated to the
fianal gff file and the path with the individual gff files creates a final
gff with the gene names modified to be more informative
'''

import argparse
import os

def replace_last_tabs_with_space(line):
    # Split the line from the right based on "\t" delimiter
    parts = line.rsplit("\t", 3)
    
    # Join the parts with a space
    updated_line = " ".join(parts)
    
    return updated_line

def gene2list(gene_file: str)-> list:
    '''
    Function to read a file with gene ids and convert it to a list

    Inputs:
        gene_file (str): path to the gene list file
    
    Outputs:
        l_genes (list): list of genes
    '''
    l_genes = list()
    with open(gene_file, "r") as f_in:
        for line in f_in:
            l_genes.append(line.strip())
    return l_genes


def concatenate(l_genes, path):
    '''
    Function that reads a list of gene id, open their corresponding gff,
    replace the name of the gene and transcript with the orthoDB id and 
    write all the gff to the same file.
    '''
    anot_fild = 2
    # Open the output file
    f_out = open("complete.gff", "w")
    # Read the genes in the list
    for gen in l_genes:
        # Write the gff of the gene
        with open(path + "/" + gen + ".gff", "r") as f_in:
            for line in f_in:
                line_l = line.split()
                # Replace the gene for the orthoDB id
                if line_l[anot_fild] == "gene":
                    gen2replace = line_l[-1]
                    line = line.replace(gen2replace, gen)
                    f_out.write(line)
                # Replace the transcript for the orthoDB id
                elif line_l[anot_fild] == "transcript":
                    transcript2replace = line_l[-1]
                    line = line.replace(transcript2replace, gen + ".t1")
                    f_out.write(line)
                else:
                    line = line.replace(gen2replace, gen)
                    line = replace_last_tabs_with_space(line)
                    f_out.write(line)
    f_out.close()

           
def main():
    # Get the list of genes
    l_genes = gene2list(snakemake.input[0])
    # Concatenate the gff files and produce the final gff
    gff_path = os.path.join(snakemake.input.busco_path,f"{snakemake.params.lineage}",
                        "augustus_output","gff")
    out_path = snakemake.output[0]
    concatenate(l_genes, gff_path,out_path)

if __name__=="__main__":
    main()