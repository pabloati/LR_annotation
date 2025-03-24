'''
Script to generate gubsets of a GeneBank file of a desired size 

The biggest subset is equal to the subset of a smaller size plus some
more genes. (The smaller subsets are subsets of the bigger subsets)
'''
import random

class gene:
    '''
    Class to store the information of a gene in GeneBank format

    Attributes:
    -----------
        locus (str): gene name plus some extra information
        features (str): CDS positions
        base_count (str): number of nucleotides
        origin (str): sequence
    
    Methods:
    --------
        __str__: convert the gene object to a str
        get_gene_id: get the gene id
    '''
    def __init__(self, locus, features, base_count, origin) -> None:
        self.locus = locus
        self.features = features
        self.base_count = base_count
        self.origin = origin


    def __str__(self) -> str:
        return self.locus + self.features + self.base_count + self.origin


    def get_gene_id(self)-> str:
        gene_line = self.features.split("\n")[-2]
        gene_id = gene_line.lstrip(' /gene="')
        gene_id = gene_id.rstrip('.t1"')
        return gene_id

def parse_gb(f_in: str)-> list:
    '''
    Function to parse a GeneBank file and convert every gene in an object
    of gene class. Save all the genes in a list
    
    Inputs:
        f_in (str): path to GeneBank file
    
    Outputs:
        l_gb (list): list of gene objects
    '''
    # List of the fild of a gene in GeneBank format
    l_campos = []
    # List of genes
    l_gb = []
    # Index of the fild that is being read
    indice = 0
    with open(f_in, "r") as gb_file:
        for linea in gb_file:
            # LOCUS is the first fild. If a list of filds already exists
            # convert the list to gene object
            if linea.startswith("LOCUS"):
                if l_campos:    
                    gene_aux = gene(l_campos[0], 
                                    l_campos[1], 
                                    l_campos[2], 
                                    l_campos[3])
                    l_gb.append(gene_aux)
                    l_campos = list()
                    indice = 0

                if not l_campos:
                    l_campos =["","","",""]
                    # add LOCUS to its corresponding position in the list
                    l_campos[indice] += linea
            # When FEATURES is read change index to 1
            elif linea.startswith("FEATURES"):
                indice = 1
                # Add FEATURES to it corresponding position in the list
                l_campos[indice] += linea
            # When BASE is read change index to 2
            elif linea.startswith("BASE"):
                indice = 2
                # Add BASE to it corresponding position in the list
                l_campos[indice] += linea
            # When ORIGIN is read change index to 3
            elif linea.startswith("ORIGIN"):
                indice = 3
                # Add ORIGIN to it corresponding position in the list
                l_campos[indice] += linea
            # When the line in the GeneBank file doesn't start with any of
            # the previous str add the line to the saved index 
            else:
                l_campos[indice] += linea

        # Save the last gene
        gene_aux = gene(l_campos[0], l_campos[1], l_campos[2], l_campos[3])
        l_gb.append(gene_aux)

        return l_gb

def main():
    # Read the GeneBank file and generate a list of gene objects
    l_gb = parse_gb(snakemake.input[0])

    # Generate a random list with the number of genes of each subset
    n_genes = min(int(snakemake.params.size), len(l_gb)) 
    
    # Set the seed
    random.seed(snakemake.params.seed)
    
    # Generate a random list of int to use as index of the list of genes
    # and sample those genes
    randomlist = random.sample(range(0, n_genes), n_genes)
    
    # Generate all the subsets
    with open(snakemake.output[0], "w") as f_out:
        for j in range(n_genes):
            # Use the random lsit to sample the genes
            f_out.write(str(l_gb[randomlist[j]]))
            
if __name__=="__main__":
    main()