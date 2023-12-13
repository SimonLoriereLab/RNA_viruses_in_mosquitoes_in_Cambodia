library(tidyverse)
library(dendextend)
library(RColorBrewer)
library(viridis)



##### function
order_pairwise_dist_from_df <- function(ska_df, 
                                        all_samples, 
                                        distance_variable, 
                                        s1, s2){
  pairwise_distances_ord_correct <- c()
  count <-  1
  count_na <- 1
  for (n in all_samples[1:(length(all_samples)-1)]){
    print(n)
    print(count)
    count <-  count + 1
    for (m in all_samples[count:length(all_samples)]){
      #print(m)
      
      pairwise_dist <- as.numeric(ska_df[which(ska_df[,s1] == n & ska_df[,s2] == m), distance_variable])
      if (is.na(pairwise_dist)){
        pairwise_dist <- as.numeric(ska_df[which(ska_df[,s2] == n & ska_df[,s1] == m),distance_variable])
        #print("-")
      }else{
        #print("+")
      }
      
      if(is.na(pairwise_dist)){
        print("NA distance")
        temp_df <- tibble(n, m)
        if (count_na  == 1){
          nadf <- temp_df
        }else{
          nadf <- bind_rows(nadf, temp_df)
        }
        count_na = count_na + 1
      }
      pairwise_distances_ord_correct <- c(pairwise_distances_ord_correct, pairwise_dist)
      #print(pairwise_distances_ord_correct)
      
    }
    
  }
  
  if (count_na  == 1){
    print("all good, no pairs with NA distance")
  }else{
    write.table(file = "/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.na_pairs____shouldnt_exist____.txt",
                x = nadf, append = F,quote = F,sep = "\t", row.names = F, col.names = F)
  }
  
  return(pairwise_distances_ord_correct)
}









ska_df <- read_delim("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.distances.tsv",
                      delim = "\t")

mdf <- read_delim("/full_path_to/wd/mosquito_species_verification/metadata/CMVM_full_metadata_updated_20211021.tsv",
                     delim = "\t")

# mdf <- read_delim("/full_path_to/wd/mosquito_species_verification/metadata/CMVM_full_metadata_updated_20211021.tsv",
#                   delim = "\t")

tdf <- read_delim("/full_path_to/wd/QC/analysis/Trimmomatic/trimmomatic_summary_CMVM_20211122.txt",
                  delim = "\t")

tdf1 <- tdf[which(tdf$prop_qualtrim_kept_PE_reads > 0.9),] ### I checked the threshold

list_of_poorly_sequenced_samples_to_remove <- c("F1073", "F1288", "F1234")

ska_df1 <- ska_df[which(!str_detect(ska_df$`Sample 1`, "PBS") &
                          !str_detect(ska_df$`Sample 2`, "PBS") &
                          !ska_df$`Sample 2` %in% list_of_poorly_sequenced_samples_to_remove &
                          !ska_df$`Sample 1` %in% list_of_poorly_sequenced_samples_to_remove), ]

nrow(ska_df)
nrow(ska_df1)



# ### THIS HAS BEEEN DONE AND PAIRS with NA SNP distances were noted and removed from the final list


# fs <- list.files("/full_path_to/wd/mosquito_species_verification/analysis/SKA/split_kmer_files/")
# countf=1
# for(f in fs){
#   print(f)
#   fp <- paste0("/full_path_to/wd/mosquito_species_verification/analysis/SKA/split_kmer_files/",f)
#   print(file.info(fp)$size)
# 
#   finfo_temp_df <- tibble(f, fsize = file.info(fp)$size)
#   if (countf==1){
#     finfo_df <- finfo_temp_df
#   }else{
#     finfo_df <- bind_rows(finfo_df, finfo_temp_df)
#   }
#   countf= countf + 1
# }
# 
# finfo_df$sample <- str_remove(finfo_df$f, "\\.skf")
# 
# finfo_df1 <- finfo_df %>% arrange(fsize)
# 




# all_samples <- sort(union(sort(unique(ska_df1$`Sample 1`)), sort(unique(ska_df1$`Sample 2`))))
# all_samples2 <- finfo_df1[which(!str_detect(finfo_df1$sample, "PBS") & 
#                                   !finfo_df1$sample %in% list_of_poorly_sequenced_samples_to_remove),]$sample
# 
# print("all_samples==all_samples2")
# sort(all_samples)==sort(all_samples2)
# sum(sort(all_samples)==sort(all_samples2))
# length(all_samples)==length(all_samples2)


# write.table(file = "/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samplenames_ordered.txt",
#             x = all_samples2, append = F,quote = F,sep = "\t", row.names = F, col.names = F)
all_samples2 <- read_lines("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samplenames_ordered.txt")
#all_samples <- read_lines("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samplenames_ordered.txt")

# ord_distances <-  order_pairwise_dist_from_df(ska_df = ska_df1, all_samples = all_samples2, distance_variable = "SNP distance", s1 = "Sample 1", s2 = "Sample 2")
# write.table(file = "/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_pairdist_ordered.txt",
#             x = ord_distances, append = F,quote = F,sep = "\t", row.names = F, col.names = F)

ord_distances <- read_lines("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_pairdist_ordered.txt")

# #ord_distances <- read_lines("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_pairdist_ordered.txt")
# 
# 
mydist = structure(ord_distances,
                   Size = length(all_samples2), Labels = all_samples2,
                   Diag = FALSE, Upper = FALSE,
                   class = "dist")


hclust(mydist) %>%
  as.dendrogram() -> dend
# 
# 
# 

# morph_sp <- sort(unique(mdf[which(mdf$Sample %in% all_samples2),]$genus_species))
# num_cols <- length(morph_sp)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_palette <- sample(col_vector, num_cols)


# 
# # 
# # plot(dend)
# # 
# # # Chart (left)
# # dend %>% 
# #   # Custom branches
# #   set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
# #   # Custom labels
# #   set("labels_col", "orange") %>% set("labels_cex", 1) %>%
# #   plot()
# # # Middle
# # dend %>% 
# #   set("nodes_pch", 19)  %>% 
# #   set("nodes_cex", 0.7) %>% 
# #   set("nodes_col", "orange") %>% 
# #   plot()
# # # right
# # dend %>% 
# #   set("leaves_pch", 19)  %>% 
# #   set("leaves_cex", 0.7) %>% 
# #   set("leaves_col", "skyblue") %>% 
# #   plot()
# # 
# # dend %>% 
# #   set("leaves_pch", 19)  %>% 
# #   set("leaves_cex", 0.7) %>% 
# #   set("leaves_col", c(1:15)) %>% 
# #   plot()
# # 
# # Color in function of the cluster
# 
# png("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_dendrogram.png",
#     width = 10000, height = 1000)
# 
# par(mar=c(1,1,1,7))
# dend %>%
#   set("labels_col", value = col_palette, k=num_cols) %>%
#   set("branches_k_color", value = col_palette, k=num_cols) %>%
#   set("nodes_cex", 2) %>%
#   set("leaves_cex", 4) %>%
#   plot(horiz=FALSE, axes=FALSE)
# 
# dev.off()





# pdf("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_dendrogram.pdf",
#     width = 80, height = 12)
# 
# par(mar=c(7,7,7,7))
# dend %>%
#   set("labels_col", value = col_palette, k=num_cols) %>%
#   set("branches_k_color", value = col_palette, k=num_cols) %>%
#   set("leaves_cex", 0.5) %>%
#   plot(horiz=FALSE, axes=FALSE)
# 
# dev.off()
# 
# 
# morph_genus <- sort(unique(mdf[which(mdf$Sample %in% all_samples2),]$genus))
# num_cols2 <- length(morph_genus)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector2 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_palette2 <- sample(col_vector2, num_cols2)
# 
# pdf("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_dendrogram.genus.pdf",
#     width = 12, height = 5)
# 
# par(mar=c(2,2,2,2))
# dend %>%
#   #set("labels_col", value = col_palette2, k=num_cols2) %>%
#   set("labels_cex", 0.05) %>%
#   #set("leaves_col", value = col_palette2, k=num_cols2) %>%
#   set("branches_k_color", value = col_palette2, k=num_cols2) %>%
#   #set("leaves_cex", 0) %>%
#   #highlight_branches_lwd %>%
#   set("branches_lwd", 1) %>%
#   plot(horiz=FALSE, axes=FALSE)
# 
# dev.off()


dendord_df <- tibble(Sample = labels(dend))
dendord_df <- left_join(dendord_df, mdf, by = "Sample") %>% select(sample_id=Sample, genus, genus_species)
dendord_df <- left_join(dendord_df, tdf1, by = "sample_id")  %>% select(sample_id, genus, genus_species, num_qualtrim_kept_PE_reads)
dendord_df$lg_kept_PE_reads <- log10(dendord_df$num_qualtrim_kept_PE_reads)

#first bar
nColor=4
col_palette_magma4 <- magma(n=nColor, direction = -1)

min(dendord_df$lg_kept_PE_reads)
max(dendord_df$lg_kept_PE_reads)

rank <- cut(dendord_df$lg_kept_PE_reads, breaks=c(4,5,6,7,8), 
                                     labels=c("4 <= 5", "5 <= 6", "6 <= 7" , "7 <= 8"), 
                                     ordered_result = T, include.lowest = TRUE)

lg_read_num <- col_palette_magma4[ rank ]


the_bars <- as.data.frame(tibble(lg_read_num))

#species bars
for (sp in sort(unique(dendord_df$genus_species))){
  temp_v <- ifelse(dendord_df$genus_species == sp, "#000000", "#f0f0f0")
  the_bars <- cbind(the_bars, temp_v)
  names(the_bars)[names(the_bars) == 'temp_v'] <- str_replace_all(sp, " ", "_")
}
the_bars=the_bars[,order(ncol(the_bars):1)]






pdf("/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_dendrogram.custom_color.1.pdf",
    width = 10, height = 7)

par(mar=c(14,3.5,0.05,0.05))
dend %>%
  #set("labels_col", value = col_palette2, k=num_cols2) %>%
  set("labels_cex", 0.1) %>%
  #set("leaves_col", value = col_palette2, k=num_cols2) %>%
  #set("branches_k_color", value = col_palette2, k=num_cols2) %>%
  #set("leaves_cex", 0) %>%
  #highlight_branches_lwd %>%
  #set("branches_lwd", 1) %>%
  plot(horiz=FALSE, axes=FALSE)
colored_bars(colors = the_bars, dend = dend, sort_by_labels_order = FALSE, cex.rowLabels = 0.5, y_scale = 0.6, y_shift = -0.04)

dev.off()



class(dend)

library(phylogram)

write.dendrogram(x = dend, 
                 file = "/full_path_to/wd/mosquito_species_verification/analysis/SKA/distances/Culicidae_transcriptome_comparison_all.all_samples_dendrogram.nwk", 
                 edges = T)
