# Script to filter the SQANTI3 output classificaiton file using information
# blast

# Load libraries
library(dplyr)
library(ggplot2)

# Args
class_path = snakemake.input.class
class_out_path = snakemake.output.class  # Output for the filtered classification
gtf_path = snakemake.input.gtf
gtf_out = snakemake.output.gtf

# Read the classification and the blast output
classification = read.delim(class_path)

# Genetal filtering of the classification file
filtered_clasification <- classification %>% 
  filter(exons > 1) %>% # Remove monoexons
  filter(coding == "coding") %>% # Keep only the coding genes
  filter(predicted_NMD == "FALSE") %>%  # Remove transcript with NMD signals
  filter(perc_A_downstream_TTS < 60) %>% # Remove intrapriming candidates
  filter(all_canonical == "canonical") %>% # Keep trancript with canonical SJ
  filter(CDS_length > 300) # Keep trancripts with a CDS of at least 300 nt

# Filtering using the short reads coverage of the SJ infroamtion   
filtered_clasification <- filtered_clasification %>%  
  filter(RTS_stage == "FALSE") %>% # Not RTS
  filter((ORF_length + 1)*3 == CDS_length) # Complete CDS with Stop codon

# To generate the hints no more filtering steps are needed, write the filtered
# classification and exit

  write.table(filtered_clasification, file = class_out_path,
            row.names = F, quote = F, sep = "\t")


gtf_columns = c("chr","source","feature","start","end","score","strand","frame","attributes")
main_gtf <- read.delim(gtf_path,header=F, col.names = gtf_columns)
main_gtf %>% rowwise() %>%
  mutate(transcript_id = gsub("transcript_id ","", strsplit(attributes,";")[[1]][1])) %>%
  filter(transcript_id %in% filtered_clasification$isoform) %>%
  select(-transcript_id) %>%
  write.table(file = gtf_out, row.names = F, sep = "\t",header = F)


