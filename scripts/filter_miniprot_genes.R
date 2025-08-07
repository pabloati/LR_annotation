library(dplyr)
library(readr)
library(stringr)

# Args
args <- commandArgs(trailingOnly = TRUE)
gff_path <- args[1]
output <- args[2]

gff_cols <- c("chr","source","type","start","end","score","strand","phase","attributes")

df <- read_tsv(gff_path, col_names = gff_cols)

df %>% mutate(id = gsub(".*=","",str_split(attributes,";", simplify = TRUE)[,1])) %>%
    group_by(id) %>% 
    summarise(type = paste0(unique(type),collapse=";"),
        attributes = first(attributes)) %>%
    ungroup() %>%
    mutate(identity = str_extract(attributes, "Identity=([^;]+)") %>% gsub("Identity=","",.) %>% as.numeric()) %>%
    filter(identity >= 0.95 & str_detect(type,"stop_codon")) %>% pull(id) -> good_genes

df %>% mutate(id = gsub(".*=","",str_split(attributes,";", simplify = TRUE)[,1])) %>%
    filter(id %in% good_genes) %>%
    select(-id) %>%
    write_tsv(output,col_names=F)
