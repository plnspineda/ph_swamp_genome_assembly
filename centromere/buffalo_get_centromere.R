library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("Usage: Rscript get_centromere.R [input.sat-ctr.repbed.out] [out name]")
}

# Repeat masker out file as input needs to be cleaned first. Machine readable with only the satellite/ctr repeats.

# setwd("/home/polen/Documents/Tuli_and_Wagyu/WxC_initialanalysis/Charolais/centromere/")
# pathfile <- "Charolais.haplotype1.chrNames.20231003_tmp_asm.sat-ctr.out.tsv"
# inname <- "Charolais"

pathfile <- args[1]
inname <- args[2]

infile <- read_tsv(pathfile, col_names = TRUE)
infile <- infile %>% unite(new_column, 15:ncol(infile), sep = "-")

names(infile) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                   "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                   "repeat_end","repeat_left","ID", "remarks", "others")

head(infile)

chr <- as.character(1:29)
chrX <- "X"
chrY <- "Y"
chr1_X <- c(chr,chrX,chrY)

infile$query_sequence[!infile$query_sequence %in% chr1_X] <- "Unplaced"

#loop thro and add query_seq_align_len and perc_id
infile_subset <- infile %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

infile_subset$assembly <- rep(inname,nrow(infile_subset))

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",infile$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",infile$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

infile$new_repeat_begin <- new_repeat_begin
infile$new_repeat_end <- infile$repeat_end
infile$new_repeat_left <- new_repeat_left

#top align length of repeat family with highest presence
infile_subset_topRepeatLen <- infile_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

### use l-apply to all the autosomes

chrom_list <- c(1:29,"X","Y")
group_func <- function(chrom) {
  df <- infile %>%
    filter(family == "Satellite/centr" & query_sequence == chrom) %>%
    mutate(query_align_len = query_end - query_begin + 1,
           perc_id = 1 - (divergence / 100),
           diff = query_end - lag(query_begin)) %>%
    filter(perc_id > 0.6) %>%
    select("query_sequence", "query_begin", "query_end", "query_left", "strand", "repeat", "family", "query_align_len", "diff")
  
  df$group <- cumsum(!is.na(df$diff) & df$diff > 1000000)
  
  df_summary <- df %>%
    group_by(group) %>%
    summarise(chromosome = chrom,
              start = first(query_begin),
              end = last(query_end),
              rep_size = sum(query_align_len),
              loc_size = last(query_end) - first(query_begin) + 1,
              asm = inname)
  
  assign(paste("HS_ctrloc_", chrom, sep = ""), df_summary, envir = .GlobalEnv)
}

invisible(lapply(chrom_list, group_func))

infile_ctr <- do.call(rbind, lapply(chrom_list, function(chrom) get(paste("HS_ctrloc_", chrom, sep = ""))))
infile_ctr_final <- infile_ctr %>%
  group_by(chromosome) %>%
  top_n(1, rep_size)
infile_plot <- ggplot(infile_ctr_final, aes(x = factor(chromosome, levels = c(1:29, "X", "Y")), y = rep_size)) +
  geom_col(position = "identity", alpha = 0.3) +
  labs(title = "Size of centromeres per chromosomes", x = "Chromosome", y = "Size (bp)")
print(infile_plot)

## plot the centromeres
infile_plot <- ggplot(infile_ctr_final, aes(x = factor(chromosome, levels = c(1:29, "X", "Y")))) +
  geom_col(aes(y = rep_size, fill = "Satellite repeats"), position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_col(aes(y = loc_size, fill = "Centromere length"), position = position_dodge(width = 0.8), alpha = 0.7) +
  labs(title = "Size of centromeres per chromosomes", x = "Chromosome", y = "Values") +
  scale_fill_manual(name = "Category", values = c("Satellite repeats" = "blue", "Centromere length" = "gray")) +
  theme(plot.title = element_text(hjust = 0.5))
png(paste0(inname,"_ctr_count.png"), width = 10, height = 8, units = "in", res = 300)
print(infile_plot)
dev.off()

ctr_type_plot <- infile_subset %>% 
  ggplot(aes(x = factor(query_sequence, levels = c(1:29, "X", "Y")), y = query_align_len, fill = `repeat`)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(title = "Centromere Length per Chromosomes", x = "Chromosomes", y = "Size", fill = "Repeat Families") +
  theme(plot.title = element_text(hjust = 0.5))
png(paste0(inname,"_ctr_reptype.png"), width = 10, height = 6, units = "in", res = 300)
print(ctr_type_plot)
dev.off()

repeat_classes <- infile_subset %>% select(`repeat`) %>% group_by(`repeat`) %>% summarise(count=n())
write_tsv(repeat_classes, file = paste0(inname,"_ctr_repeatclasses.txt"))

infile_ctr_final <- infile_ctr_final %>% mutate(perc = rep_size/loc_size*100)
write_tsv(infile_ctr_final, file = paste0(inname,"_ctr_allchr.txt"))

infile_ctr_final_filter <- infile_ctr_final %>% filter(rep_size > 100000, perc > 50)
write_tsv(infile_ctr_final_filter, file = paste0(inname,"_ctr_final.txt"))
