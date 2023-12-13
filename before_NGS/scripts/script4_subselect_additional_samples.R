getwd()
setwd("/full_path_to/wd/before_NGS/")

library(tidyverse)
library(readxl)
library(viridis)
library(measurements)
library(ggrepel)
library(gridExtra)
library(cowplot)

all_df <- read.delim2("metadata/CMVM - priority annotated full sample list 26072021.tsv",
                      sep = "\t", colClasses = "character")
all_df$Batch <- as.numeric(all_df$Batch)
done_df <- all_df[which(all_df$Batch > 0),]
yet_df <- all_df[which(all_df$Batch == 0),]

unique(yet_df$genus_species)
### SELECT Anopheles separately
Anopheles_df <- yet_df[which(yet_df$genus_species == "Anopheles peditaeniatus"),]
yet_df1 <- yet_df[which(!yet_df$genus_species == "Anopheles peditaeniatus"),]

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
#yet_df1_sub <- subsample_n_fun(df = yet_df1, traits = c('genus_species', 'Location'), n = 4)
#yet_df1_sub_and_Anoph <- bind_rows(yet_df1_sub, Anopheles_df)
# write.table(x = yet_df1_sub_and_Anoph, file = "metadata/subsampled_for_suplementary_sequencing.txt",
#            append = F, quote = F, sep = "\t", col.names = T, row.names = F)

##### Manually subselect further
yet_df1_sub_final <- read.delim2("metadata/subsampled_for_suplementary_sequencing.txt",
                               sep = "\t", colClasses = "character") %>% 
  arrange(Location, new_box, new_eppendorf)

yet_df1_sub_final$Batch <- "19"
yet_df1_sub_final$LibPool <- "21"
yet_df1_sub_final[31:60,'Batch'] <- "20"
yet_df1_sub_final[31:60,'LibPool'] <- "22"

yet_df$Batch <- as.character(yet_df$Batch)

yet_df2 <- bind_rows(yet_df1_sub_final, yet_df[which(!yet_df$Sample %in% yet_df1_sub_final$Sample),])

done_df$Batch <- as.character(done_df$Batch)

all_df1 <- bind_rows(done_df, yet_df2)

# write.table(x = all_df1, file = "priority_annotated_full_sample_list_20210905.txt", 
#             append = F, quote = F, sep = "\t'", col.names = T, row.names = F)


