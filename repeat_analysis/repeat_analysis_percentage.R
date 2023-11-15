library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)

# tail -n +3 RM_WB.out
# grep -v '^$' RM_WB.out > tmp
# sed -e 's/  */\t/g' tmp > tmp1
# sed 's/^\t//' < tmp1 > tmp2
# mv tmp2 RM_WB.out

dir1 <- "/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/swamp_buffalo/repeat_analysis/Rscript/data/"
figures <- "/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/swamp_buffalo/repeat_analysis/Rscript/figures/"
output <- "/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/swamp_buffalo/repeat_analysis/Rscript/output/"

#read in philippine swamp buffalo
path1 <- paste0(dir1,"RM_SB_withunplaced.out")
pcc_uoa_sb_1 <- read_tsv(path1, col_names = FALSE)

names(pcc_uoa_sb_1) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                         "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                         "repeat_end","repeat_left","ID", "remarks")

chr <- as.character(1:24)
chrX <- "X"
chr1_X <- c(chr,chrX)

pcc_uoa_sb_1$query_sequence[!pcc_uoa_sb_1$query_sequence %in% chr1_X] <- "Unplaced"

#loop thro and add query_seq_align_len and perc_id
pcc_uoa_sb_1_subset <- pcc_uoa_sb_1 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

pcc_uoa_sb_1_subset$assembly <- rep("PCC_UOA_SB_1v2",nrow(pcc_uoa_sb_1_subset))

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",pcc_uoa_sb_1$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",pcc_uoa_sb_1$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

pcc_uoa_sb_1$new_repeat_begin <- new_repeat_begin
pcc_uoa_sb_1$new_repeat_end <- pcc_uoa_sb_1$repeat_end
pcc_uoa_sb_1$new_repeat_left <- new_repeat_left

#top align length of repeat family with highest presence
pcc_uoa_sb_1_subset_topRepeatLen <- pcc_uoa_sb_1_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

################################################################################

#read in italian mediterranean buffalo
path2 <- paste0(dir1,"RM_WB_withunplaced.out")
uoa_wb_1 <- read_tsv(path2, col_names = FALSE)

names(uoa_wb_1) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                     "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                     "repeat_end","repeat_left","ID", "remarks")

uoa_wb_1$query_sequence[!uoa_wb_1$query_sequence %in% chr1_X] <- "Unplaced"

#loop thro and add query_seq_align_len and perc_id
uoa_wb_1_subset <- uoa_wb_1 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

uoa_wb_1_subset$assembly <- rep("UOA_WB_1",nrow(uoa_wb_1_subset))

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",uoa_wb_1$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",uoa_wb_1$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

uoa_wb_1$new_repeat_begin <- new_repeat_begin
uoa_wb_1$new_repeat_end <- uoa_wb_1$repeat_end
uoa_wb_1$new_repeat_left <- new_repeat_left

#top align length of repeat family with highest presence
uoa_wb_1_subset_topRepeatLen <- uoa_wb_1_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

##################################################################################

#read in indian murrah buffalo
path3 <- paste0(dir1,"NDDB_SH_1.out")
nddb_sh_1 <- read_tsv(path3, col_names = FALSE)

names(nddb_sh_1) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                      "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                      "repeat_end","repeat_left","ID", "remarks")

nddb_sh_1$query_sequence[!nddb_sh_1$query_sequence %in% chr1_X] <- "Unplaced"

#loop thro and add query_seq_align_len and perc_id
nddb_sh_1_subset <- nddb_sh_1 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

nddb_sh_1_subset$assembly <- rep("NDDB_SH_1",nrow(nddb_sh_1_subset))

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",nddb_sh_1$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",nddb_sh_1$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

nddb_sh_1$new_repeat_begin <- new_repeat_begin
nddb_sh_1$new_repeat_end <- nddb_sh_1$repeat_end
nddb_sh_1$new_repeat_left <- new_repeat_left

#top align length of repeat family with highest presence
nddb_sh_1_subset_topRepeatLen <- nddb_sh_1_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

#########################################################

#read in chinese swamp buffalo
path4 <- paste0(dir1,"CUSA_SWP_withunplaced.out")
cusa_swp <- read_tsv(path4, col_names = FALSE)

names(cusa_swp) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                     "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                     "repeat_end","repeat_left","ID", "remarks")

cusa_swp$query_sequence[!cusa_swp$query_sequence %in% chr1_X] <- "Unplaced"

#loop thro and add query_seq_align_len and perc_id
cusa_swp_subset <- cusa_swp %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

cusa_swp_subset$assembly <- rep("CUSA_SWP",nrow(cusa_swp_subset))

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",cusa_swp$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",cusa_swp$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

cusa_swp$new_repeat_begin <- new_repeat_begin
cusa_swp$new_repeat_end <- cusa_swp$repeat_end
cusa_swp$new_repeat_left <- new_repeat_left

#top align length of repeat family with highest presence
cusa_swp_subset_topRepeatLen <- cusa_swp_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

##################################################################################

#read in chinese river buffalo
path5 <- paste0(dir1,"CUSA_RVB_frozen_renamed.fna.out")
cusa_rvb <- read_tsv(path5, col_names = FALSE)

names(cusa_rvb) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                     "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                     "repeat_end","repeat_left","ID", "remarks")

cusa_rvb$query_sequence[!cusa_rvb$query_sequence %in% chr1_X] <- "Unplaced"

#loop thro and add query_seq_align_len and perc_id
cusa_rvb_subset <- cusa_rvb %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

cusa_rvb_subset$assembly <- rep("CUSA_RVB",nrow(cusa_rvb_subset))

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",cusa_rvb$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",cusa_rvb$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

cusa_rvb$new_repeat_begin <- new_repeat_begin
cusa_rvb$new_repeat_end <- cusa_rvb$repeat_end
cusa_rvb$new_repeat_left <- new_repeat_left

#top align length of repeat family with highest presence
cusa_rvb_subset_topRepeatLen <- cusa_rvb_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

#read in male swamp buffalo
path6 <- paste0(dir1,"Male_swamp_chr_CRA007045.renamed.fna.out")
male_swamp <- read_tsv(path6, col_names = FALSE)

names(male_swamp) <- c("score", "divergence", "deletion", "insertion", "query_sequence", "query_begin", 
                       "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                       "repeat_end","repeat_left","ID", "remarks")

chr <- as.character(1:25)
chrX <- "X"
chrY <- "Y"
chr1_XY <- c(chr,chrX,chrY)

male_swamp$query_sequence[!male_swamp$query_sequence %in% chr1_X] <- "Unplaced"

#loop thro and add query_seq_align_len and perc_id
male_swamp_subset <- male_swamp %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

male_swamp_subset$assembly <- rep("Wang_2023",nrow(male_swamp_subset))

#repeat_begin and repeat_left have silly ()
new_repeat_begin <- gsub("\\(","",male_swamp$repeat_begin)
new_repeat_begin <- gsub("\\)","",new_repeat_begin)
new_repeat_begin <- as.integer(new_repeat_begin)

new_repeat_left <- gsub("\\(","",male_swamp$repeat_left)
new_repeat_left <- gsub("\\)","",new_repeat_left)
new_repeat_left <- as.integer(new_repeat_left)

male_swamp$new_repeat_begin <- new_repeat_begin
male_swamp$new_repeat_end <- male_swamp$repeat_end
male_swamp$new_repeat_left <- new_repeat_left

#top align length of repeat family with highest presence
male_swamp_subset_topRepeatLen <- male_swamp_subset %>% 
  group_by(family) %>% summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))

################################################################################

##################################################################################

#Combine assembly datasets
pcc_uoa_sb_1_subset_noUnplaced <- pcc_uoa_sb_1_subset %>% filter(query_sequence != "Unplaced")
uoa_wb_1_subset_noUnplaced <- uoa_wb_1_subset %>% filter(query_sequence != "Unplaced")
nddb_sh_1_subset_noUnplaced <- nddb_sh_1_subset %>% filter(query_sequence != "Unplaced")
cusa_swp_subset_noUnplaced <- cusa_swp_subset %>% filter(query_sequence != "Unplaced")
cusa_rvb_subset_noUnplaced <- cusa_rvb_subset %>% filter(query_sequence != "Unplaced")
male_swamp_subset_noUnplaced <- male_swamp_subset %>% filter(query_sequence != "Unplaced")

all_spp_rm_subset <- rbind(pcc_uoa_sb_1_subset_noUnplaced, male_swamp_subset_noUnplaced, uoa_wb_1_subset_noUnplaced, nddb_sh_1_subset_noUnplaced, cusa_swp_subset_noUnplaced, cusa_rvb_subset_noUnplaced)


#order spp
order_spp <- c("PCC_UOA_SB_1v2", "Wang_2023", "UOA_WB_1", "NDDB_SH_1", "CUSA_SWP", "CUSA_RVB")
all_spp_rm_subset$assembly <- factor(all_spp_rm_subset$assembly, levels = order_spp)

all_spp_rm_subset_nofil_summ <- all_spp_rm_subset %>% group_by(family,assembly) %>% summarise(count = n()) %>%
  arrange(desc(count))

all_spp_rm_subset_summ <- all_spp_rm_subset %>% group_by(family,assembly) %>% summarise(count = n()) %>% 
  arrange(desc(count)) %>% filter(family == "LINE/L1" | family == "LINE/RTE-BovB" | family == "Satellite/centr")

#### trying violin plot to show query_align_len distribution 
all_spp_rm_subset_violin <- all_spp_rm_subset %>% group_by(assembly) %>% 
  filter(family == "LINE/L1" | family == "LINE/RTE-BovB" | family == "Satellite/centr") %>%
  filter(query_align_len > 2.5e3)

#The shape of violin starts and ends with the IQR ends
png(filename = paste0(figures,"RepeatsFamiliyByLengthWaterBuffalo.png"), width = 12, height = 8, units = "in", res = 300)
bp <- ggplot(all_spp_rm_subset_violin, aes(x=assembly, y=query_align_len, group=assembly)) 
bp <- bp + geom_violin(aes(fill=assembly)) + guides(fill=FALSE) + geom_boxplot(width = 0.1, outlier.shape = NA) ## add median and quartile
bp <- bp + scale_y_log10(labels = scales::comma,breaks=c(0,5000,10000,15000,20000,25000,30000,40000,50000,60000,80000,100000,120000,150000,250000,500000))
bp <- bp + facet_grid(. ~ family) + theme_bw(base_size = 17) 
bp <- bp + theme(strip.text.x = element_text(size = 20), axis.text.x = element_text(angle = 70, hjust = 1), element_line(colour = "black"))
bp <- bp + ylab("Length of reference matched to repeat (bp)") + xlab("Assembly")
bp
dev.off()

####################
#for supplementary

#Combine pcc_uoa_sb_1 and uoa_wb_1 subsets to look at unplaced vs in chr
all_spp_rm_subset_supp <- rbind(pcc_uoa_sb_1_subset, uoa_wb_1_subset, nddb_sh_1_subset, cusa_swp_subset, cusa_rvb_subset, male_swamp_subset)

all_spp_rm_subset_supp_summ <- all_spp_rm_subset_supp %>% group_by(query_sequence,family,assembly) %>% summarise(count = n()) %>%
  filter(family == "LINE/L1" | family == "LINE/RTE-BovB" | family == "Satellite/centr") %>%
  mutate(placement = ifelse(query_sequence == "Unplaced","Unplaced","Chromosome")) %>% arrange(desc(count))

#barplot of different repeat family count by sequence placement status
png(filename =  paste0(figures,"RepeatsFamilyswamp_river.png"),width = 12, height = 4, units = "in", res = 300)
h <- ggplot(data = all_spp_rm_subset_supp_summ, aes(x = reorder(as.factor(family), -count), y = log10(count), fill = assembly)) + 
  geom_bar(stat = "identity", position="dodge") + ylab(expression(log[10]*" count")) + xlab("Repeat family")
h <- h + facet_grid(~placement) 
h <- h + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 0, hjust = 0.5))
h <- h + theme_bw(base_size = 12)
h
dev.off()

all_spp_rm_subset_violin_summary <- all_spp_rm_subset_violin %>% group_by(assembly, `family`) %>% summarise(query_align_len_sum = sum(query_align_len))

##### <START> PERCENTAGE OF REPEATS in the GENOME #####
#loop swamp and calc align query length but dont filter percent id to get %repeat
pcc_uoa_sb_1_subset_nofilpercid <- pcc_uoa_sb_1 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

#swamp repeat of all classes straight from repeat masker/genome size
perc_pcc <- sum(pcc_uoa_sb_1_subset_nofilpercid$query_align_len)/2899768601
#0.5145273

#####

#loop Mediterranean and calc align query length but dont filter percent id to get %repeat
uoa_wb_1_subset_nofilpercid <- uoa_wb_1 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

#river repeat of all classes straight from repeat masker/genome size
perc_uoa <- sum(uoa_wb_1_subset_nofilpercid$query_align_len)/2655776128
#0.4841492 (paper is 47.48%)

#####

#loop indian murrah and calc align query length but dont filter percent id to get %repeat
nddb_sh_1_subset_nofilpercid <- nddb_sh_1 %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

#river repeat of all classes straight from repeat masker/genome size
perc_nddb <- sum(nddb_sh_1_subset_nofilpercid$query_align_len)/2622460639
#0.4823201

#####

#loop chinese swamp and calc align query length but dont filter percent id to get %repeat
cusa_swp_subset_nofilpercid <- cusa_swp %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")

#swamp repeat of all classes straight from repeat masker/genome size
perc_cswap <- sum(cusa_swp_subset_nofilpercid$query_align_len)/2631349264
#0.4805031

#loop chinese swamp and calc align query length but dont filter percent id to get %repeat
cusa_rvb_subset_nofilpercid <- cusa_rvb %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")
#swamp repeat of all classes straight from repeat masker/genome size
perc_crvb <- sum(cusa_rvb_subset_nofilpercid$query_align_len)/2645519677
#0.4805031

#loop male swamp and calc align query length but dont filter percent id to get %repeat
male_swamp_subset_nofilpercid <- male_swamp %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  select("query_sequence","query_begin","query_end","query_left","strand","repeat","family","query_align_len")
#swamp repeat of all classes straight from repeat masker/genome size
perc_maleswamp <- sum(male_swamp_subset_nofilpercid$query_align_len)/2899768601


##### <END> PERCENTAGE OF REPEATS in the GENOME #####

##### <START> Highest repeat family in the unplaced #####
pcc_uoa_sb_1_subset_Unplaced <- pcc_uoa_sb_1_subset %>% filter(query_sequence == "Unplaced") 
pcc_uoa_sb_1_subset_Unplaced_topRepeatLen <- pcc_uoa_sb_1_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
perc_unplace_pcc <- pcc_uoa_sb_1_subset_Unplaced_topRepeatLen$total_len[1]/sum(pcc_uoa_sb_1_subset_Unplaced_topRepeatLen$total_len)
# 131243073 bases of Satellite/centr in the unplaced, which is the highest in terms of length.
# 0.8349815 of the repeats in the unplaced are Satellite/centr

uoa_wb_1_subset_Unplaced <- uoa_wb_1_subset %>% filter(query_sequence == "Unplaced")
uoa_wb_1_subset_Unplaced_topRepeatLen <- uoa_wb_1_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
perc_unplace_uoa <- uoa_wb_1_subset_Unplaced_topRepeatLen$total_len[1]/sum(uoa_wb_1_subset_Unplaced_topRepeatLen$total_len)
# 131243073 bases of Satellite/centr in the unplaced, which is the highest in terms of length.
# 0.5422177 of the repeats in the unplaced are Satellite/centr

nddb_sh_1_subset_Unplaced <- nddb_sh_1_subset %>% filter(query_sequence == "Unplaced")
nddb_sh_1_subset_Unplaced_topRepeatLen <- nddb_sh_1_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
perc_unplace_nddb <- nddb_sh_1_subset_Unplaced_topRepeatLen$total_len[1]/sum(nddb_sh_1_subset_Unplaced_topRepeatLen$total_len)
#  bases of Satellite/centr in the unplaced, which is the highest in terms of length.
#  of the repeats in the unplaced are Satellite/centr

cusa_swp_subset_Unplaced <- cusa_swp_subset %>% filter(query_sequence == "Unplaced")
cusa_swp_subset_Unplaced_topRepeatLen <- cusa_swp_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
perc_unplace_cswp <- cusa_swp_subset_Unplaced_topRepeatLen$total_len[1]/sum(cusa_swp_subset_Unplaced_topRepeatLen$total_len)
# 15482451 bases of Satellite/centr in the unplaced, which is the highest in terms of length.
# 0.3701227 of the repeats in the unplaced are Satellite/centr 

cusa_rvb_subset_Unplaced <- cusa_rvb_subset %>% filter(query_sequence == "Unplaced")
cusa_rvb_subset_Unplaced_topRepeatLen <- cusa_rvb_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
perc_unplace_crvb <- cusa_rvb_subset_Unplaced_topRepeatLen$total_len[1]/sum(cusa_rvb_subset_Unplaced_topRepeatLen$total_len)
#  bases of Satellite/centr in the unplaced, which is the highest in terms of length.
#  of the repeats in the unplaced are Satellite/centr 

male_swamp_subset_Unplaced <- male_swamp_subset %>% filter(query_sequence == "Unplaced") 
male_swamp_subset_Unplaced_topRepeatLen <- male_swamp_subset_Unplaced %>% group_by(family) %>% 
  summarise(total_len = sum(query_align_len)) %>% 
  arrange(desc(total_len))
perc_unplace_maleswamp <- male_swamp_subset_Unplaced_topRepeatLen$total_len[1]/sum(male_swamp_subset_Unplaced_topRepeatLen$total_len)
# 131243073 bases of Satellite/centr in the unplaced, which is the highest in terms of length.
# 0.8349815 of the repeats in the unplaced are Satellite/centr

### saving files to the folder

perc_repeat_genomes <- c("PCC_UOA_SB_1v2", "Wang_2023", "CUSA_SWP", "UOA_WB_1", "NDDB_SH_1", "CUSA_RVB")
perc_repeat_value <- c(perc_pcc, perc_maleswamp, perc_cswap, perc_uoa, perc_nddb, perc_crvb)
perc_repeat_unplaced <- c(perc_unplace_pcc, perc_unplace_maleswamp, perc_unplace_cswp, perc_unplace_uoa, perc_unplace_nddb,  perc_unplace_crvb )
perc_repeat_df <- data.frame(perc_repeat_genomes, perc_repeat_value, perc_repeat_unplaced)

write_tsv(pcc_uoa_sb_1_subset_topRepeatLen, file = paste0(output,"pcc_uoa_sb_1_subset_topRepeatLen.txt"))
write_tsv(uoa_wb_1_subset_topRepeatLen, file = paste0(output,"uoa_wb_1_subset_topRepeatLen.txt"))
write_tsv(nddb_sh_1_subset_topRepeatLen, file = paste0(output,"nddb_sh_1_subset_topRepeatLen.txt"))
write_tsv(cusa_swp_subset_topRepeatLen, file = paste0(output,"cusa_swp_subset_topRepeatLen.txt"))
write_tsv(cusa_rvb_subset_topRepeatLen, file = paste0(output,"cusa_rvb_subset_topRepeatLen.txt"))
write_tsv(male_swamp_subset_topRepeatLen, file = paste0(output,"male_swamp_subset_topRepeatLen.txt"))

write_tsv(pcc_uoa_sb_1_subset_Unplaced_topRepeatLen, file = paste0(output,"pcc_uoa_sb_1_subset_Unplaced_topRepeatLen.txt"))
write_tsv(uoa_wb_1_subset_Unplaced_topRepeatLen, file = paste0(output,"uoa_wb_1_subset_Unplaced_topRepeatLen.txt"))
write_tsv(nddb_sh_1_subset_Unplaced_topRepeatLen, file = paste0(output,"nddb_sh_1_subset_Unplaced_topRepeatLen.txt"))
write_tsv(cusa_swp_subset_Unplaced_topRepeatLen, file = paste0(output,"cusa_swp_subset_Unplaced_topRepeatLen.txt"))
write_tsv(cusa_rvb_subset_Unplaced_topRepeatLen, file = paste0(output,"cusa_rvb_subset_Unplaced_topRepeatLen.txt"))
write_tsv(male_swamp_subset_Unplaced_topRepeatLen, file = paste0(output,"male_swamp_subset_Unplaced_topRepeatLen.txt"))

write_tsv(perc_repeat_df, file = paste0(output,"RepeatPercentageGenomes.txt"))
write_tsv(all_spp_rm_subset_violin, file = paste0(output,"all_spp_rm_subset_violin.txt"))
write_tsv(all_spp_rm_subset_violin_summary, file = paste0(output,"all_spp_rm_subset_violin_summary.txt"))

print("Finished.")