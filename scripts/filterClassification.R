# Script to filter the SQANTI3 output classificaiton file using information
# blast

# Load libraries
library(dplyr)
library(ggplot2)

# Args
args = commandArgs(trailingOnly=TRUE)
class_path = snakemake.input[0]
class_out_path = snakemake.output[0]  # Output for the filtered classification
blast_file = args[4] # Blast output
tecnologia = args[5] # Sequencing tech

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
  filter((ORF_length + 1)*3 == CDS_length) # COmplete CDS with Stop codon

# To generate the hints no more filtering steps are needed, write the filtered
# classification and exit
if (length(args) != 5){
  write.table(filtered_clasification, file = class_out_path,
            row.names = F, quote = F, sep = "\t")
  quit()
}

# Set thresholds based on tech
th <- 0.75 # if the technology is ONT or MIX use the 3th quertile

# if the tech is PB use the median
if (tecnologia == "PB"){
  th <- 0.5
}

# Filtering based on the SJ coverage
filtered_clasification <- filtered_clasification %>%  
  filter(min_cov > quantile(filtered_clasification$min_cov, th)) # Minimun coverage of the SJ

# Read the blast output
blast <- read.delim(blast_file, header = F)

# Prepare the blast output to be used
blast$qcoverage <- blast$V13/blast$V14 # Query coverage
blast$seq_hit <- paste(blast$V1, blast$V2, sep = "_") # Id of the query hit duo
sel <- match(unique(blast$seq_hit), blast$seq_hit) # Keep only one query hit duo
blast <- blast[sel,]
# Summary the stats of the 3 best hits of the query
blast_top3h_summary <- blast %>% 
  group_by(V1) %>% # group by query id
  summarise(n_hit = n(), # number of hits
            qcover = max(qcoverage), # Select the maximun query coverage of the 3 hits
            identidad = max(V3), min(V11)) # select the maximun identity %
# Reduce the number of hits in case of having more than 3
blast_top3h_summary$n_hit <- ifelse(blast_top3h_summary$n_hit > 3, 
                                    3, blast_top3h_summary$n_hit)

# Plot the query coverage considering the number of hits
blast_plot <- ggplot(blast_top3h_summary, 
                     aes(x=qcover, group=as.factor(n_hit), 
                     fill=as.factor(n_hit))) +
  geom_density(alpha = 0.5) + 
  geom_vline(xintercept =  0.85, linetype = "longdash") + 
  xlim(c(0,1.5)) +
  guides(fill=guide_legend(title="Nº of blastp hits")) + 
  xlab("Query coverage")

# Keep those blast hits with a query coverage higher than and 85% and lower 
# than a 120%
blast_85 <- blast_top3h_summary %>% 
  filter(qcover > 0.85) %>% 
  filter(qcover < 1.2)

# if the technology is PB the classification is filteres using only the blast
# hits. If the technology is ONT then the FL counts are also used. Initial 
# test showed that filtering only with blast using PB that was enough to get
# high values of precision. To get comparable levels of precision using 
# nanopore it is also needed to filter using FL counts
if (tecnologia == "PB"){
  filtered_clasification <- filtered_clasification %>% 
    filter(isoform %in% blast_85$V1)
} else{
  filtered_clasification <- filtered_clasification %>% 
    filter(isoform %in% blast_85$V1) %>% 
    filter(FL > quantile(classification$FL, th))
}

# Plot the number of exons per transcript before and after
p_exon_pre <- ggplot(classification, aes(exons))+geom_histogram(binwidth=1) + 
  geom_histogram(binwidth=1,fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  xlab("Number of exons") + ylab("Number of transcripts") + geom_vline(xintercept = median(classification$exons))
p_exon_post <-  ggplot(filtered_clasification, aes(exons))+geom_histogram(binwidth=1) + 
  geom_histogram(binwidth=1,fill="#69b3a2", color="#e9ecef", alpha=0.9) + 
  xlab("Number of exons") + ylab("Number of transcripts") + geom_vline(xintercept = median(filtered_clasification$exons))

# Write the output
class_out_path = paste(out_dir, paste0("filtered_",sq_prefix, "_classification.txt"), sep = "/")
write.table(filtered_clasification, file = class_out_path,
            row.names = F, quote = F, sep = "\t")

ggsave(paste(out_dir,"exon_pre.png", sep = "/"), p_exon_pre, width = 7, height = 5)
ggsave(paste(out_dir,"exon_post.png", sep = "/"), p_exon_post, width = 7, height = 5)
ggsave(paste(out_dir,"blast.png", sep = "/"), blast_plot, width = 7, height = 5)