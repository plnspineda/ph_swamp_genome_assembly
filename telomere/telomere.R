rm(list = ls())
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggpubr)
# library(plyr)
library(RColorBrewer)

## read data

# SB_size <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/fai_files/PCC_UOA_SB_1_frozen.fna.fai", col_names = FALSE)
# names(SB_size) <- c("id","size","start_loc","A","B")
# SB_size <- SB_size %>% select(id, size) %>% mutate(endloc = cumsum(size), begin = endloc - size +1, stop = endloc)
# SB2_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/PCC_UOA_SB_1_chronly.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)
# # SB_telorev_subset <- SB_telomere %>% group_by(id) %>% top_n(1, reverse_repeat_number)
# # SB_telofor_subset <- SB_telomere %>% group_by(id) %>% top_n(1, forward_repeat_number)

SB_size <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/fai_files/PCC_UOA_SB_1_frozen.fna.fai", col_names = FALSE)
names(SB_size) <- c("id","size","start_loc","A","B")
SB_size <- SB_size %>% select(id, size) %>% mutate(size = size + 4000000, endloc = cumsum(size), begin = endloc - size)
SB2_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/PCC_UOA_SB_1_chronly.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)

WB_size <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/fai_files/UOA_WB_1_frozen.fna.fai", col_names = FALSE)
names(WB_size) <- c("id","size","start_loc","A","B")
WB_size <- WB_size %>% select(id, size) %>% mutate(size = size + 4000000, endloc = cumsum(size), begin = endloc - size)
WB_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/UOA_WB_1_chronly.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)

ND_size <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/fai_files/NDDB_SH_1_frozen.fna.fai", col_names = FALSE)
names(ND_size) <- c("id","size","start_loc","A","B")
ND_size <- ND_size %>% select(id, size) %>% mutate(size = size + 4000000, endloc = cumsum(size), begin = endloc - size)
ND_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/NDDB_SH_1_frozen.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)

CHS_size <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/fai_files/CUSA_SWP_chronly_renamed.fna.fai", col_names = FALSE)
names(CHS_size) <- c("id","size","start_loc","A","B")
CHS_size <- CHS_size %>% select(id, size) %>% mutate(size = size + 4000000, endloc = cumsum(size), begin = endloc - size)
CHS_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/CUSA_SWP_chronly_renamed.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)

CHR_size <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/fai_files/CUSA_RVB_chronly_renamed.fna.fai", col_names = FALSE)
names(CHR_size) <- c("id","size","start_loc","A","B")
CHR_size <- CHR_size %>% select(id, size) %>% mutate(size = size + 4000000, endloc = cumsum(size), begin = endloc - size)
CHR_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/CUSA_RVB_chronly_renamed.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)

MLE_size <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/fai_files/Male_swamp_CRA007045_chronly_renamed.fna.fai", col_names = FALSE)
names(MLE_size) <- c("id","size","start_loc","A","B")
MLE_size <- MLE_size %>% select(id, size) %>% mutate(size = size + 4000000, endloc = cumsum(size), begin = endloc - size)
MLE_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/Male_swamp_chr_CRA007045.renamed.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)


# ARS_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/other_species/ARS-UCD1.3_chronly.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)
# T2TWa_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/other_species/Wagyu_sire_cattle.haplotype2.fasta_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)
# T2THu_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/other_species/T2T-CHM13v2.0_frozen.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)
# Gor_telomere <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/other_species/NHGRI_mGorGor1-v1-0.2_chronly.fna_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)



## filter the non-telomere sites
SB_telo_filtered <- SB2_telomere %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>%
  filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X")

test <- SB_telo_filtered %>% group_by(id, tmp = window == 20000) %>% 
  summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
  ungroup() %>% select(-tmp)

SB_telo_filtered <- merge(SB_telo_filtered, SB_size, by="id") %>% 
  mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
         loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
         telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
         telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
         telomere_size = telomere_count*6, asm="PCC_UOA_SB_1")
SB_telo_final <- SB_telo_filtered %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)

WB_telo_filtered <- WB_telomere %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>% 
  filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X") %>% group_by(id, tmp = window == 20000) %>% 
  summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
  ungroup() %>% select(-tmp)
WB_telo_filtered <- merge(WB_telo_filtered, WB_size, by="id") %>% 
  mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
         loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
         telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
         telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
         telomere_size = telomere_count*6, asm="UOA_WB_1")
WB_telo_final <- WB_telo_filtered %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)

ND_telo_filtered <- ND_telomere %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>% 
  filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X") %>% group_by(id, tmp = window == 20000) %>% 
  summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
  ungroup() %>% select(-tmp)
ND_telo_filtered <- merge(ND_telo_filtered, ND_size, by="id") %>% 
  mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
         loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
         telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
         telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
         telomere_size = telomere_count*6, asm="NDDB_SH_1")
ND_telo_final <- ND_telo_filtered %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)

CHS_telo_filtered <- CHS_telomere %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>% 
  filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X") %>% group_by(id, tmp = window == 20000) %>% 
  summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
  ungroup() %>% select(-tmp)
CHS_telo_filtered <- merge(CHS_telo_filtered, CHS_size, by="id") %>% 
  mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
         loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
         telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
         telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
         telomere_size = telomere_count*6, asm="CUSA_SWP")
CHS_telo_final <- CHS_telo_filtered %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)

CHR_telo_filtered <- CHR_telomere %>% group_by(id) %>%slice(c(1:2, (n()-1):n())) %>% 
  filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X") %>% group_by(id, tmp = window == 20000) %>% 
  summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
  ungroup() %>% select(-tmp)
CHR_telo_filtered <- merge(CHR_telo_filtered, CHR_size, by="id") %>% 
  mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
         loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
         telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
         telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
         telomere_size = telomere_count*6, asm="CUSA_RVB")
CHR_telo_final <- CHR_telo_filtered %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)

MLE_telo_filtered <- MLE_telomere %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>% 
  filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X", id!="Y") %>% group_by(id, tmp = window == 20000) %>% 
  summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
  ungroup() %>% select(-tmp)
MLE_telo_filtered <- merge(MLE_telo_filtered, MLE_size, by="id") %>% 
  mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
         loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
         telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
         telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
         telomere_size = telomere_count*6, asm="Wang_2023")
MLE_telo_final <- MLE_telo_filtered %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)

## list sizes 
all_telo_size <- genome_list <- list(PCC_UOA_SB_1 = SB_telo_final$telomere_size,
                                CUSA_SWP_1 = CHS_telo_final$telomere_size,
                                UOA_WB_1 = WB_telo_final$telomere_size,
                                CUSA_RVB = CHR_telo_final$telomere_size,
                                NDDB_SH_1 = ND_telo_final$telomere_size,
                                Wang_2023 = MLE_telo_final$telomere_size)

all_asm <- rbind(SB_telo_final, CHS_telo_final, WB_telo_final, CHR_telo_final, ND_telo_final, MLE_telo_final)

## for supplementary table

PCC_telo_sup <- SB_telo_filtered  %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
  arrange(as.numeric(id))
UOA_telo_sup <- WB_telo_filtered %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
  arrange(as.numeric(id))
NDD_telo_sup <- ND_telo_filtered %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
  arrange(as.numeric(id))
CHR_telo_sup <- CHR_telo_filtered %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
  arrange(as.numeric(id))
CHS_telo_sup <- CHS_telo_filtered %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
  arrange(as.numeric(id))
MLE_telo_sup <- MLE_telo_filtered %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
  arrange(as.numeric(id))
all_telo_sup <- rbind(PCC_telo_sup, UOA_telo_sup, NDD_telo_sup, CHR_telo_sup, CHS_telo_sup, MLE_telo_sup)

chr <- data.frame(c(1:24))
names(chr) <- "id"
all_asm$id <- as.numeric(all_asm$id)
all_asm_chr <- merge(chr, all_asm, by="id", all = TRUE) %>% select(id, telomere_count, asm) %>%
  group_by(id, asm) %>% summarise(total_telomere_count = sum(telomere_count))

## for supplementary table

telomere_counts_per_chr <- all_asm_chr %>% spread(key = asm, value = total_telomere_count, fill = NA) %>% select(id, PCC_UOA_SB_1, Wang_2023, NDDB_SH_1, UOA_WB_1, CUSA_SWP, CUSA_RVB)
names(telomere_counts_per_chr) <- c("Chromosomes", "PCC_UOA_SB_1v2", "Wang_2023", "NDDB_SH_1", "UOA_WB_1", "CUSA_SWP", "CUSA_RVB")
write_tsv(telomere_counts_per_chr, file = "/home/polen/Documents/swamp_buffalo/centromere/tidk/TelomereCountsPerAssemblies.-withmaletxt")

#telomere mean

all_asm_compute_mean <- all_asm_chr %>% na.omit() %>% group_by(asm) %>% summarise_at(vars(total_telomere_count),list(mean = mean, sd = sd)) %>%
  as.data.frame()
telomere_mean_sd <- ggplot(all_asm_compute_mean, aes(x=asm, y=mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
  geom_point() +
  labs(title = "Telomere size in water buffalo assemblies", x = "Genome assembies", y = "Mean (unit)") +
  theme(plot.title = element_text(hjust = 0.5))
png("/home/polen/Documents/swamp_buffalo/centromere/tidk/telomere_mean_sd_withmale.png", width = 10, height = 8, units = "in", res = 300)
print(telomere_mean_sd)
dev.off()

### for bargraph

asm <- c("PCC_UOA_SB_1v2", "NDDB_SH_1", "UOA_WB_1", "CUSA_SWP", "CUSA_RVB", "Wang_2023")
colors <- c("PCC_UOA_SB_1v2"="#ff6361", "NDDB_SH_1"="#27aeef","UOA_WB_1"="#00bfa0","CUSA_SWP"="#ffa600","CUSA_RVB"="#bc5090", "Wang_2023"="#6C3483")
all_asm_size_per_chromosome <- telomere_counts_per_chr %>% gather(key = "Assemblies", value = "Unit", -Chromosomes)
all_asm_size_per_chromosome_plot <- all_asm_size_per_chromosome %>%
  ggplot(aes(x = factor(Chromosomes), y = Unit, fill = Assemblies)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of Telomere Counts", x = "Chromosomes", y = "Unit") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = colors, breaks = asm)
png("/home/polen/Documents/swamp_buffalo/centromere/tidk/all_asm_size_per_chromosome_plot.png", width = 10, height = 6, units = "in", res = 300)
print(all_asm_size_per_chromosome_plot)
dev.off()

### notes
# 1. to make a gap between chromosomes in the bedgraph
## a. identify the last line of each group
## b. add another line: start with the last value of the previous number + 40000 with 0 value

## plot telomeres per chromosomes
gray_colors <- rep(c("#A9A9A9", "#D3D3D3"), length.out = 24)


### this works
# SB_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/PCC_UOA_SB_1_chronly.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
# names(SB_bedgraph) <- c("id","start","end","value")
# SB_bedgraph <- SB_bedgraph %>% filter(id != "X")
# SB_chromsize <- SB_size %>% head(n=23)
# SB_bedgraph_2 <- left_join(SB_bedgraph, SB_chromsize, by = "id")
# SB_bedgraph_final <- SB_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
# SB_bedgraph_final$id <- factor(SB_bedgraph_final$id, levels = c(1:23))
# SB_bedgraph_rects <- SB_chromsize %>% select(id, begin, stop) %>% group_by(id) %>% summarise(begin = min(begin), stop = max(stop)) %>% arrange(as.numeric(id))
# SB_bedgraph_rects$id <- factor(SB_bedgraph_rects$id, levels = c(1:23))
# SB_bg_plot <- ggplot() +
#   geom_rect(data = SB_bedgraph_rects, aes(xmin = begin, xmax = stop, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
#   geom_step(data = SB_bedgraph_final, aes(x = start, y = value, color = id)) +
#   scale_fill_manual(values = gray_colors) +
#   theme_classic() +
#   ylim(0, 2000) +
#   xlim(0,2600000000) +
#   xlab("") +
#   ylab("")
# print(SB_bg_plot)


########################### test ###########################

# SB_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/PCC_UOA_SB_1_chronly.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
# names(SB_bedgraph) <- c("id","start","end","value")
# SB_bedgraph <- SB_bedgraph %>% filter(id != "X")
# 
# # to add new row at the end of the group
# SB_bedgraph <- SB_bedgraph %>%
#   group_by(id) %>%
#   group_modify(~ .x %>%
#                  summarise(start = max(start) + 40000, end = max(end) + 40000, value = 0) %>%
#                  add_row(.x, .)
#   ) %>%
#   ungroup()
# 
# SB_chromsize <- SB_size %>% head(n=23) #%>% mutate(begin = begin+40000)
# SB_bedgraph_2 <- left_join(SB_bedgraph, SB_size, by = "id")
# SB_bedgraph_final <- SB_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
# SB_bedgraph_final$id <- factor(SB_bedgraph_final$id, levels = c(1:23))
# SB_bedgraph_rects <- SB_chromsize %>% select(id, begin, stop) %>% group_by(id) %>% summarise(begin = min(begin), stop = max(stop)) %>% arrange(as.numeric(id))
# SB_bedgraph_rects$id <- factor(SB_bedgraph_rects$id, levels = c(1:23))
# SB_bg_plot <- ggplot() +
#   geom_rect(data = SB_bedgraph_rects, aes(xmin = begin, xmax = stop, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
#   geom_step(data = SB_bedgraph_final, aes(x = start, y = value, color = id)) +
#   scale_fill_manual(values = gray_colors) +
#   theme_classic() +
#   ylim(0, 2000) +
#   xlim(0,2600000000) +
#   xlab("") +
#   ylab("")
# print(SB_bg_plot)

#################################################################################


## trying to add gaps between chromosomes
# test <- SB_bedgraph %>% group_by(id) %>% slice(n()) %>% mutate(start = last(end), end = last(end)+200000, value = 0) %>% ungroup
# test_out <- bind_rows(SB_bedgraph, test) %>% arrange(id, start)
# SB_chromsize <- SB_size %>% head(n=23)
# SB_bedgraph_2 <- left_join(test_out, SB_chromsize, by = "id")
# SB_bedgraph_final <- SB_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
# SB_bedgraph_final$id <- factor(SB_bedgraph_final$id, levels = c(1:23))
# SB_bedgraph_rects <- SB_chromsize %>% select(id, begin, stop) %>% group_by(id) %>% summarise(begin = min(begin), stop = max(stop)) %>% arrange(as.numeric(id))
# SB_bedgraph_rects$id <- factor(SB_bedgraph_rects$id, levels = c(1:23))
# SB_bg_plot <- ggplot() +
#   geom_rect(data = SB_bedgraph_rects, aes(xmin = begin, xmax = stop, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
#   geom_step(data = SB_bedgraph_final, aes(x = start, y = value, color = id)) +
#   scale_fill_manual(values = gray_colors) +
#   theme_classic() +
#   ylim(0, 2000) +
#   xlim(0,2600000000) +
#   xlab("") +
#   ylab("") +
#   print(SB_bg_plot)

## trying to add gaps between chromosomes
# SB_chromsize <- SB_size %>% head(n=23)
# SB_bedgraph_2 <- left_join(SB_bedgraph, SB_chromsize, by = "id")
# SB_bedgraph_final <- SB_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
# SB_bedgraph_final$id <- factor(SB_bedgraph_final$id, levels = c(1:23))
# SB_bedgraph_rects <- SB_chromsize %>% select(id, begin, stop) %>% group_by(id) %>% summarise(begin = min(begin), stop = max(stop)) %>% arrange(as.numeric(id))
# SB_bedgraph_rects$id <- factor(SB_bedgraph_rects$id, levels = c(1:23))
# SB_bg_plot <- ggplot() +
#   geom_rect(data = SB_bedgraph_rects, aes(xmin = begin, xmax = stop, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
#   geom_step(data = SB_bedgraph_final, aes(x = start, y = value, color = id)) +
#   scale_fill_manual(values = gray_colors) +
#   theme_classic() +
#   ylim(0, 2000) +
#   xlim(0,2600000000) +
#   xlab("") +
#   ylab("") +
# print(SB_bg_plot)
# SB_bg_plot <-ggplot(SB_bedgraph_final, aes(x=start, y=value, color=id)) +
#   geom_step() +
#   theme_classic() +
#   ylim(0, 2000) + xlab("") + ylab("")
# print(SB_bg_plot)

SB_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/PCC_UOA_SB_1_chronly.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
names(SB_bedgraph) <- c("id","start","end","value")
SB_bedgraph <- SB_bedgraph %>% filter(id != "X")

# to add new row at the end of the group
SB_bedgraph <- SB_bedgraph %>%
  group_by(id) %>%
  group_modify(~ .x %>%
                 summarise(start = max(start) + 2000000, end = max(end) + 2000000, value = 0) %>%
                 add_row(.x, .)
  ) %>%
  ungroup()

SB_chromsize <- SB_size %>% head(n=23) #%>% mutate(begin = begin+40000)
SB_bedgraph_2 <- left_join(SB_bedgraph, SB_size, by = "id")
SB_bedgraph_final <- SB_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
SB_bedgraph_final$id <- factor(SB_bedgraph_final$id, levels = c(1:23))
SB_bedgraph_rects <- SB_chromsize %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
SB_bedgraph_rects$id <- factor(SB_bedgraph_rects$id, levels = c(1:23))
SB_bg_plot <- ggplot() +
  geom_rect(data = SB_bedgraph_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
  geom_step(data = SB_bedgraph_final, aes(x = start, y = value, color = id)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2600000000) +
  xlab("") +
  ylab("")
print(SB_bg_plot)

WB_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/UOA_WB_1_chronly.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
names(WB_bedgraph) <- c("id","start","end","value")
WB_bedgraph <- WB_bedgraph %>% filter(id != "X")

# to add new row at the end of the group
WB_bedgraph <- WB_bedgraph %>%
  group_by(id) %>%
  group_modify(~ .x %>%
                 summarise(start = max(start) + 2000000, end = max(end) + 2000000, value = 0) %>%
                 add_row(.x, .)
  ) %>%
  ungroup()


WB_chromsize <- WB_size %>% head(n=24)
WB_bedgraph_2 <- left_join(WB_bedgraph, WB_chromsize, by = "id")
WB_bedgraph_final <- WB_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
WB_bedgraph_final$id <- factor(WB_bedgraph_final$id, levels = c(1:24))
names(WB_bedgraph_final) <- c("Chromosome","start","end","value")
WB_bedgraph_rects <- WB_chromsize %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
WB_bedgraph_rects$id <- factor(WB_bedgraph_rects$id, levels = c(1:24))
names(WB_bedgraph_rects) <- c("Chromosome","begin","stop")
WB_bg_plot <- ggplot() +
  geom_rect(data = WB_bedgraph_rects, aes(xmin = begin, xmax = stop, ymin = -Inf, ymax = Inf, fill = Chromosome), alpha = 0.4) +
  geom_step(data = WB_bedgraph_final, aes(x = start, y = value, color = Chromosome)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2600000000) +
  xlab("") +
  ylab("")
print(WB_bg_plot)

# WB_bedgraph_chr1 <- WB_bedgraph_2 %>% filter(id=="1") %>% mutate(start_new = start+start_loc, end_new = end+start_loc) %>% select(id,start=start_new, end=end_new,value)
# WB_bg_chr1_plot <-ggplot(WB_bedgraph_chr1, aes(x=start, y=value)) +
#   geom_step() +
#   theme_classic() +
#   labs(title = "Telomeric signal of Chromosome 1 of UOA_WB_1", x = "Region", y = "Telomeric count") +
#   theme(plot.title = element_text(hjust = 0.5))
# png("/home/polen/Documents/swamp_buffalo/centromere/tidk/WB_bg_chr1_plot.png", width = 8, height = 3, units = "in", res = 300)
# print(WB_bg_chr1_plot)
# dev.off()

ND_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/NDDB_SH_1_frozen.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
names(ND_bedgraph) <- c("id","start","end","value")
ND_bedgraph <- ND_bedgraph %>% filter(id != "X")

# to add new row at the end of the group
ND_bedgraph <- ND_bedgraph %>%
  group_by(id) %>%
  group_modify(~ .x %>%
                 summarise(start = max(start) + 2000000, end = max(end) + 2000000, value = 0) %>%
                 add_row(.x, .)
  ) %>%
  ungroup()


ND_chromsize <- ND_size %>% head(n=24)
ND_bedgraph_2 <- left_join(ND_bedgraph, ND_chromsize, by = "id")
ND_bedgraph_final <- ND_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
ND_bedgraph_final$id <- factor(ND_bedgraph_final$id, levels = c(1:24))
ND_bedgraph_rects <- ND_chromsize %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
ND_bedgraph_rects$id <- factor(ND_bedgraph_rects$id, levels = c(1:24))
ND_bg_plot <- ggplot() +
  geom_rect(data = ND_bedgraph_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
  geom_step(data = ND_bedgraph_final, aes(x = start, y = value, color = id)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2600000000) +
  xlab("") +
  ylab("")
print(ND_bg_plot)

CHS_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/CUSA_SWP_chronly_renamed.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
names(CHS_bedgraph) <- c("id","start","end","value")
CHS_chromsize <- CHS_size %>% head(n=23)
CHS_bedgraph <- CHS_bedgraph %>% filter(id != "X")
CHS_bedgraph_2 <- left_join(CHS_bedgraph, CHS_chromsize, by = "id")
CHS_bedgraph_final <- CHS_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
CHS_bedgraph_final$id <- factor(CHS_bedgraph_final$id, levels = c(1:23))
CHS_bedgraph_rects <- CHS_chromsize %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
CHS_bedgraph_rects$id <- factor(CHS_bedgraph_rects$id, levels = c(1:23))
CHS_bg_plot <- ggplot() +
  geom_rect(data = CHS_bedgraph_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
  geom_step(data = CHS_bedgraph_final, aes(x = start, y = value, color = id)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2600000000) +
  xlab("") +
  ylab("")
print(CHS_bg_plot)

CHR_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/CUSA_RVB_chronly_renamed.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
names(CHR_bedgraph) <- c("id","start","end","value")
CHR_bedgraph <- CHR_bedgraph %>% filter(id != "X")
CHR_chromsize <- CHR_size %>% head(n=24)
CHR_bedgraph_2 <- left_join(CHR_bedgraph, CHR_chromsize, by = "id")
CHR_bedgraph_final <- CHR_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
CHR_bedgraph_final$id <- factor(CHR_bedgraph_final$id, levels = c(1:24))
CHR_bedgraph_rects <- CHR_chromsize %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
CHR_bedgraph_rects$id <- factor(CHR_bedgraph_rects$id, levels = c(1:24))
CHR_bg_plot <- ggplot() +
  geom_rect(data = CHR_bedgraph_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
  geom_step(data = CHR_bedgraph_final, aes(x = start, y = value, color = id)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2600000000) +
  xlab("") +
  ylab("")
print(CHR_bg_plot)

MLE_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/Male_swamp_chr_CRA007045.renamed.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
names(MLE_bedgraph) <- c("id","start","end","value")
MLE_chromsize <- MLE_size %>% head(n=23)
MLE_bedgraph <- MLE_bedgraph %>% filter(id != "X", id != "Y")
MLE_bedgraph_2 <- left_join(MLE_bedgraph, MLE_chromsize, by = "id")
MLE_bedgraph_final <- MLE_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
MLE_bedgraph_final$id <- factor(MLE_bedgraph_final$id, levels = c(1:23))
MLE_bedgraph_rects <- MLE_chromsize %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
MLE_bedgraph_rects$id <- factor(MLE_bedgraph_rects$id, levels = c(1:23))
MLE_bg_plot <- ggplot() +
  geom_rect(data = MLE_bedgraph_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
  geom_step(data = MLE_bedgraph_final, aes(x = start, y = value, color = id)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2600000000) +
  xlab("") +
  ylab("")
print(MLE_bg_plot)

telomere_bg_plot <- ggarrange(SB_bg_plot, MLE_bg_plot, ND_bg_plot, WB_bg_plot, CHS_bg_plot, CHR_bg_plot, 
          labels = c("PCC_UOA_SB_v1.2", "Wang_2023", "NDDB_SH_1", "UOA_WB_1" ,"CUSA_SWP", "CUSA_RVB"), font.label = list(size = 12),
          hjust=-0.8,
          common.legend = TRUE, legend = "right",
          ncol = 1, nrow = 5)
telomere_bg_plot_final <- annotate_figure(telomere_bg_plot, top = text_grob("Telomere signals in the water buffalo assemblies", face = "bold", size = 18),
                left = "Count of TTAGGG/CCCTAA", bottom = "Region")
png("/home/polen/Documents/swamp_buffalo/centromere/tidk/telomere_bg_plot_final_withmale2.png", width = 10, height = 12, units = "in", res = 300)
print(telomere_bg_plot_final)
dev.off()

## for paper figure

SB_bg_final <- SB_bg_plot + theme(legend.position = "none") +labs(title = "PCC_UOA_SB_1v2")
MLE_bg_final <- MLE_bg_plot + theme(legend.position = "none") +labs(title = "Wang_2023")
ND_bg_final <- ND_bg_plot + theme(legend.position = "none") +labs(title = "NDDB_SH_1")
WB_bg_final <- WB_bg_plot + theme(legend.position = "none") +labs(title = "UOA_WB_1")
legend_b <- get_legend(WB_bg_final + theme(legend.position = "right"))

telomere_bg_plot_paper <- ggarrange(SB_bg_final, ND_bg_final, WB_bg_final, MLE_bg_final,
                              hjust=-0.8,
                              legend.grob = legend_b, legend = "right",
                              ncol = 1, nrow = 4)
telomere_bg_plot_final_paper <- annotate_figure(telomere_bg_plot_paper, top = text_grob("Telomeric signals in the water buffalo assemblies", face = "bold", size = 18),
                                          left = "Count of TTAGGG/CCCTAA", bottom = "Region")
png("/home/polen/Documents/swamp_buffalo/centromere/tidk/telomere_bg_plot_final_paper_withmale2.png", width = 12, height = 7, units = "in", res = 300)
print(telomere_bg_plot_final_paper)
dev.off()


### another plot
CS_bg_final <- CHS_bg_plot + theme(legend.position = "none") +labs(title = "CUSA_SWP")
CR_bg_final <- CHR_bg_plot + theme(legend.position = "none") +labs(title = "CUSA_RVB")

legend_b <- get_legend(WB_bg_final + theme(legend.position = "right"))

telomere_bg_plot_paper <- ggarrange(SB_bg_final, ND_bg_final, WB_bg_final, CS_bg_final, CR_bg_final,
                                    hjust=-0.8,
                                    legend.grob = legend_b, legend = "right",
                                    ncol = 1, nrow = 5)
telomere_bg_plot_final_paper <- annotate_figure(telomere_bg_plot_paper, top = text_grob("Telomeric signals in the water buffalo assemblies", face = "bold", size = 18),
                                                left = "Count of TTAGGG/CCCTAA", bottom = "Region")
png("/home/polen/Documents/swamp_buffalo/centromere/tidk/telomere_bg_plot_final_5wb_paper2.png", width = 12, height = 12, units = "in", res = 300)
print(telomere_bg_plot_final_paper)
dev.off()


# ### testing

# SB_bedgraph <- read_delim("/home/polen/Documents/swamp_buffalo/centromere/tidk/tidk_out_2k_window/PCC_UOA_SB_1_chronly.fna_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = FALSE)
# names(SB_bedgraph) <- c("id","start","end","value")
# SB_bedgraph <- SB_bedgraph %>% filter(id != "X")
# SB_chromsize <- SB_size %>% head(n=23)
# SB_bedgraph_2 <- left_join(SB_bedgraph, SB_chromsize, by = "id")
# SB_bedgraph_final <- SB_bedgraph_2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
# SB_bedgraph_rects <- SB_chromsize %>% select(id, begin, stop) %>% group_by(id) %>% summarise(begin = min(begin), stop = max(stop)) %>% arrange(as.numeric(id))
# 
# gray_colors <- rep(c("#A9A9A9", "#D3D3D3", "#808080", "#C0C0C0"), length.out = 23)
# ggplot() +
#   geom_rect(data = SB_bedgraph_rects, aes(xmin = begin, xmax = stop, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
#   geom_step(data = SB_bedgraph_final, aes(x = start, y = value, color = id)) +
#   scale_fill_manual(values = gray_colors) +  # Use different shades of gray for geom_rect()
#   theme_classic() +
#   ylim(0, 2000) +
#   xlab("") +
#   ylab("")

#### test for normal distributions

