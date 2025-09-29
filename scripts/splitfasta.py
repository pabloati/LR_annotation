from Bio import SeqIO
import os

def read_fasta(fastafile):
    records = list(SeqIO.parse(fastafile, "fasta"))
    return records

def split_fasta(records, outdir, type):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for record in records:
        filename = os.path.join(outdir, record.id + f".{type}")
        with open(filename, "w") as f:
            SeqIO.write(record, f, type)
        f.close()

if __name__ == "__main__":
    fasta_file = snakemake.input[0]
    outdir = os.path.dirname(snakemake.output[0])
    records = read_fasta(fasta_file)
    split_fasta(records, outdir, "fasta") 
