# Script to filter the SQANTI3 output classificaiton file using information
# Load libraries
suppressMessages(library(dplyr))
library(readr)

# Args
# class_path <- as.character(snakemake@input[[1]])
# class_out_path <-  as.character(snakemake@output[[2]])  # Output for the filtered classification
# gtf_path <- as.character(snakemake@input[[1])
# gtf_out <- as.character(snakemake@output[[2]])

args <- commandArgs(trailingOnly = TRUE)
class_path <- args[1]
gtf_path <- args[2]
class_out_path <- args[3]
gtf_out_path <- args[4]

# Read the classification and the blast output
classification <- read.delim(class_path)

# Genetal filtering of the classification file
filtered_clasification <- classification %>% 
  filter(RTS_stage == "FALSE") %>% # Not RTS
  filter(exons > 1) %>% # Remove monoexons
  filter(perc_A_downstream_TTS < 60) %>% # Remove intrapriming candidates
  filter(all_canonical == "canonical")  # Keep trancript with canonical SJ

  #filter(coding == "coding") %>% # Keep only the coding genes
  #filter(predicted_NMD == "FALSE") %>%  # Remove transcript with NMD signals
  #filter(CDS_length > 300) # Keep trancripts with a CDS of at least 300 nt
  #filter((ORF_length + 1)*3 == CDS_length) # Complete CDS with Stop codon

# To generate the hints no more filtering steps are needed, write the filtered
# classification and exit

write.table(filtered_clasification, file = class_out_path,
            row.names = F, quote = F, sep = "\t")



gtf_columns <- c("chr","source","feature","start","end","score","strand","frame","attributes")
 read_tsv(gtf_path,quote="",col_names=gtf_columns,show_col_types = F) %>% rowwise() %>% 
  mutate(transcript_id = gsub("\"","",gsub("transcript_id ","", strsplit(attributes,";")[[1]][1]))) %>%
  filter(transcript_id %in% filtered_clasification$isoform) %>%
  select(-transcript_id) %>%
  write.table(file = gtf_out_path, row.names = F,col.names = F, quote = F, sep = "\t")

