getwd()
setwd("/full_path_to/wd/before_NGS/")

library(tidyverse)
library(readxl)
library(viridis)
library(measurements)
library(ggrepel)
library(gridExtra)
library(cowplot)

all_df <- read.delim2("metadata/CMVM - priority annotated full sample list 01062021.tsv",
                      sep = "\t", colClasses = "character")

lpdf <- all_df[which(all_df$seq_priority == "low"),]
hpdf <- all_df[which(all_df$seq_priority == "high"),]

unique(lpdf$genus_species)
### remove Anopheles
lpdf1 <- lpdf[which(!lpdf$genus_species == "Anopheles peditaeniatus"),]

subsample_n_fun <- function(df, traits, n){
  #### NOTE: used traits should be of class 'character'
  # Make all trait value combinations if > 1 trait is provided
  count = 0
  for (tr in traits){
    df$tr <- df[, tr]
    if (count > 0){
      df$trait <- paste0(df$trait, "_", df$tr)
    }else{
      df$trait <- df$tr
    }
    count = count + 1
  }
  
  # Subsample according to the concatenated trait
  count2 = 0
  for (t in unique(df$trait)){
    print(t)
    if (nrow(df[which(df$trait == t),]) > n){
      tdf <- df[which(df$trait == t),] %>% group_by(trait) %>% sample_n(n)
    }else{
      tdf <- df[which(df$trait == t),]
    }
    if (count2 > 0){
      tdff <- bind_rows(tdff, tdf)
    }else{
      tdff <- tdf
    }
    count2 = count2 + 1
  }
  return(tdff)
}


### Run subsampling only once
# lpdf1_sub <- subsample_n_fun(df = lpdf1, traits = c('genus_species', 'Location'), n = 5)
# write.table(x = lpdf1_sub, file = "metadata/subsampled_low_priority_species.txt", 
#             append = F, quote = F, sep = "\t'", col.names = T, row.names = F)

##### Manually subselect further
lpdf1_sub_final <- read.delim2("metadata/CMVM_bioinformatics - low_priority_selection.tsv",
                               sep = "\t", colClasses = "character") %>% 
  arrange(Location, new_box, new_eppendorf)
lpdf1_sub_final$Batch <- "18"
lpdf1_sub_final$LibPool <- "18"

lpdf2 <- bind_rows(lpdf1_sub_final, lpdf[which(!lpdf$Sample %in% lpdf1_sub_final$Sample),])

hpdf1 <- bind_rows(hpdf, lpdf2)

# write.table(x = hpdf1, file = "metadata/priority_annotated_full_sample_list_23062021.txt", 
#             append = F, quote = F, sep = "\t'", col.names = T, row.names = F)


