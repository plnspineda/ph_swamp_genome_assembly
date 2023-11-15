### PSPINEDA 2023
### for counting the number of SVs from a nucmer alignment (delta) to assemblytics .bed file
### without removing duplicates

library(readr)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("Usage: Rscript assemblytics_SVsize.R [input.bed]")
}

print(getwd())

sv_bed <- args[1]

df <- read_tsv(sv_bed)
names(df) <- gsub(x = names(df), pattern = "#", replacement = "")
df$uniqueName <- paste(df$reference, df$ref_start, df$ref_stop, df$size,
                                   df$strand, df$type, sep="_")

table(df$method,df$type)
n_occur_df <- data.frame(table(df$uniqueName))
n_occur_df[n_occur_df$Freq >1,]
dupid_df <- df[df$uniqueName %in% n_occur_df$Var1[n_occur_df$Freq >1],]$ID
df <- df[!df$ID %in% dupid_df,]
#for paper, total size of swamp buffalo affected SV
sum(df$size)
