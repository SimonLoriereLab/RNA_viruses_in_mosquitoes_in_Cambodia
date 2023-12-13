library(tidyverse)
library(ggtree)
library(treeio)
library(seqinr)
library(lubridate)

################## Original tree
# tree_file_dirty_tips <- "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree.tre"
# ### Substitute all special characters except parentheses "()" and colons ":", commas "," and periods "." by underscores "_"
# tree_text <- read_file(tree_file_dirty_tips)
# cat(paste0(sort(unique(str_extract_all(tree_text, "[[:punct:]]")[[1]])), ""))
# tree_text1 <- str_replace_all(str_replace_all(tree_text, "[[\\-\\'\\/\\[\\]]]", "_"), " ", "_")
# 
# write_file(x = tree_text1,
#            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.tre",
#            append = F)
###

############## Reimport with corrected names
# tree_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.tre"
# tree <- read.newick(file = tree_file)
# #### Quick tree check
# ggtree(tree)
# 

#### In FigTree midpoint rooted, increasing tip order
tree_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean_in_figtree.nwk"
tree <- read.newick(file = tree_file)
#### Quick tree check
ggtree(tree)


### Parse tips
tip_labels_raw <- tree$tip.label
tips_df <- tibble(tip_labels_raw, 
                  rdrp_scan_phylum = str_extract(tip_labels_raw, "[[:alpha:]]+viricota"),
                  rdrp_scan_class = str_extract(tip_labels_raw, "[[:alpha:]]+viricetes"),
                  rdrp_scan_order = str_extract(tip_labels_raw, "[[:alpha:]]+virales"),
                  rdrp_scan_family = str_extract(tip_labels_raw, "[[:alpha:]]+viridae"))


tips_df$tip_labels <- str_remove(tips_df$tip_labels_raw, "^\\'")
tips_df$tip_labels <- str_remove(tips_df$tip_labels, "\\'$")
tips_df$tip_labels <- str_replace_all(str_replace_all(tips_df$tip_labels, "[[\\-\\'\\/]]", "_"), " ", "_")


R_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/R_orf_df_curated_with_chimeric_contig_id.txt", 
                   delim = "\t")
##### get rid of duplicates
temp <- R_df[which(str_detect(R_df$orf_id, "F1621_NODE_24_length_3179_cov_5.991997")),]
nrow(temp) == 1
R_df <- R_df[which(!(R_df$orf_id == "F1621_NODE_24_length_3179_cov_5.991997_1concat2" & R_df$evalue_be == max(temp$evalue_be))),]
temp <- R_df[which(str_detect(R_df$orf_id, "F1621_NODE_24_length_3179_cov_5.991997")),]
nrow(temp) == 1

temp <- R_df[which(str_detect(R_df$orf_id, "F1358_NODE_2309_length_1852_cov_1.580412")),]
nrow(temp) == 1
R_df <- R_df[which(!(R_df$orf_id == "F1358_NODE_2309_length_1852_cov_1.580412_1concat2" & is.na(R_df$species))),]
temp <- R_df[which(str_detect(R_df$orf_id, "F1358_NODE_2309_length_1852_cov_1.580412")),]
nrow(temp) == 1

R_df_reduced <- R_df %>% select(tip_labels = orf_id, 
                                contig_id, contig_set, 
                                target)


R_BE_tax_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/best_hit_taxonomy.txt", 
                          delim = "\t")

R_BE_tax_df_reduced <- R_BE_tax_df %>% select(contig_id,
                                              rdrp_scan_vir_group, 
                                              BestHit_phylum = phylum, 
                                              BestHit_class = class, 
                                              BestHit_order = order,
                                              BestHit_family = family,
                                              BestHit_genus = genus,
                                              BestHit_species = species,
                                              BestHit_pident = pident,
                                              suggested_vir_sp)

R_df_reduced <- left_join(R_df_reduced, R_BE_tax_df_reduced, by = "contig_id")



tips_df$status <- "reference"
tips_df[which(tips_df$tip_labels %in% R_df_reduced$tip_labels),]$status <- "novel"


tips_df1 <- left_join(tips_df, R_df_reduced)



# ggt4 <- ggtree(tree, layout="daylight")
# save(ggt4,
#      file = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_daylight.rdata")

load("/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_daylight.rdata")

tree_tip_df <- tibble(tip_labels_raw)
tree_tip_df1 <- left_join(tree_tip_df, tips_df1, by  = "tip_labels_raw")

tree_tip_df1$phylum <- tree_tip_df1$rdrp_scan_phylum
tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_phylum)),]$phylum <- tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_phylum)),]$BestHit_phylum 

############# PLOT


ggt41 <- ggt4  
ggt41$layers[[1]]$geom$default_aes$size <- 0
ggt41$layers[[1]]$geom$default_aes$alpha <- 0


ggt6 <- ggt41 %<+% tree_tip_df1

ggt6$data[which(is.na(ggt6$data$status)),]$status <- "internal_branch"
ggt6$data[which(is.na(ggt6$data$status)),]$status <- factor(ggt6$data[which(is.na(ggt6$data$status)),]$status, 
                                                            levels = c("internal_branch", "reference", "novel" ), ordered = T)



check_df <- ggt6$data


ggt6 <- ggt6 + geom_tree(mapping = aes(color = status, alpha = status, size = status), layout = "daylight")
ggt6 <- ggt6 + geom_tippoint(aes(shape = phylum, color = status), size = 1.5, alpha = 0.75)
ggt6 <- ggt6 + geom_treescale()


ggt6 <- ggt6 + scale_color_manual(values = c("internal_branch"= "black", "reference"= "gray70","novel" = "firebrick3"))
ggt6 <- ggt6 + scale_alpha_manual(values = c("internal_branch"= 0.75,"reference"= 0.75,"novel" = 0.75))
ggt6 <- ggt6 + scale_size_manual(values = c("internal_branch"= 0.5,"reference"= 0.5,"novel" = 1.5))




ggt6 <- ggt6 + theme(legend.position = "bottom", legend.direction = "vertical")
print(ggt6)

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_daylight_2.pdf", 
       plot = ggt6, device = "pdf", height = 22, width = 20, units = "cm")

#############



#############



#############



#############



#############



#############



#############












##### TEMP plot to annotate by nodes
ggt41 <- ggt4  
ggt41$layers[[1]]$geom$default_aes$size <- 0
ggt41$layers[[1]]$geom$default_aes$alpha <- 0


ggt6 <- ggt41 %<+% tree_tip_df1

ggt6$data[which(is.na(ggt6$data$status)),]$status <- "internal_branch"
ggt6$data[which(is.na(ggt6$data$status)),]$status <- factor(ggt6$data[which(is.na(ggt6$data$status)),]$status, 
                                                            levels = c("internal_branch", "reference", "novel" ), ordered = T)

ggt6$data$rdrp_scan_detected_vir_group <- NA
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Pisuviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Pisuviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Kitrinoviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Kitrinoviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Duplornaviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Duplornaviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Lenarviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Lenarviricota"

ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Monjiviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Monjiviricetes"
ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Ellioviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Ellioviricetes"
ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Insthoviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Insthoviricetes"

ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Ghabrivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Ghabrivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Sobelivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Sobelivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Durnavirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Durnavirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Amarillovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Amarillovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Wolframvirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Wolframvirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Ourlivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Ourlivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Martellivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Martellivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Picornavirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Picornavirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Reovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Reovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Tolivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Tolivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Tymovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Tymovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Nidovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Nidovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Nodamuvirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Nodamuvirales"

ggt6$data[which(ggt6$data$rdrp_scan_family == regex("Permutotetraviridae", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Permutotetraviridae"

ggt6$data[which(ggt6$data$status == "internal_branch"),]$rdrp_scan_detected_vir_group <- "internal_branch"


ggt6 <- ggt6 + geom_tree(mapping = aes(color = rdrp_scan_detected_vir_group, alpha = status, size = status), layout = "daylight")
#ggt6 <- ggt6 + geom_label2(aes(label=node), fill='lightgreen', size = 1, alpha = 1/4)
#ggt6 <- ggt6 + geom_tippoint(aes(fill = phylum), shape = 21, size = 5, alpha = 0.9)
ggt6 <- ggt6 + geom_tippoint(aes(size = rdrp_scan_detected_vir_group), shape = 21, alpha = 0.9)


ggt6 <- ggt6 + geom_treescale()


### ggt6 <- ggt6 + scale_color_manual(values = c("internal_branch"= "black", "reference"= "gray70","novel" = "firebrick3"))
ggt6 <- ggt6 + scale_alpha_manual(values = c("internal_branch"= 0.75,"reference"= 0.75,"novel" = 0.75))
ggt6 <- ggt6 + scale_size_manual(values = c("internal_branch"= 0.5,"reference"= 0.5,"novel" = 1.5))

# ggt6 <- ggt6 + scale_fill_manual(values = c("Duplornaviricota" = "purple", "Kitrinoviricota"= "green", 
#                                              "Lenarviricota"= "yellow", "Negarnaviricota"= "red", 
#                                              "Pisuviricota" = "blue"))

# ggt6 <- ggt6 + scale_fill_manual(values = c('Caliciviridae' = ,'Dicistroviridae','Durnavirales','Iflaviridae','Nidovirales','Partitiviridae','Patatavirales',
#                                             'Picornavirales','Picornaviridae','Polycipiviridae','Sobelivirales','Stellavirales',
#                                             
#                                             
#                                             
#                                             ))


# cat(paste0(sort(unique(ggt6$data[which(ggt6$data$rdrp_scan_phylum == "Pisuviricota"),]$rdrp_scan_LCA)), collapse = "','"))


ggt6 <- ggt6 + theme(legend.position = "bottom", legend.direction = "vertical")
#print(ggt6)

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_daylight_2_temp.pdf", 
       plot = ggt6, device = "pdf", height = 100, width = 100, units = "cm")














##### Another tree form
ggt41 <- ggt4  
ggt41$layers[[1]]$geom$default_aes$size <- 0
ggt41$layers[[1]]$geom$default_aes$alpha <- 0


ggt6 <- ggt41 %<+% tree_tip_df1

ggt6$data[which(is.na(ggt6$data$status)),]$status <- "internal_branch"
ggt6$data[which(is.na(ggt6$data$status)),]$status <- factor(ggt6$data[which(is.na(ggt6$data$status)),]$status, 
                                                            levels = c("internal_branch", "reference", "novel" ), ordered = T)


#### only taxa for which novel sequences were detected
ggt6$data$rdrp_scan_detected_vir_group <- NA
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Pisuviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Pisuviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Kitrinoviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Kitrinoviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Duplornaviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Duplornaviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Lenarviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Lenarviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Negarnaviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Negarnaviricota"



temp_df <- ggt6$data[which(!ggt6$data$rdrp_scan_phylum %in% c("Pisuviricota", "Kitrinoviricota", "Duplornaviricota", "Lenarviricota", "Negarnaviricota") & !ggt6$data$status == "internal_branch"),]





ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Monjiviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Monjiviricetes"
ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Ellioviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Ellioviricetes"
ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Insthoviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Insthoviricetes"

ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Ghabrivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Ghabrivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Sobelivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Sobelivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Durnavirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Durnavirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Amarillovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Amarillovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Wolframvirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Wolframvirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Ourlivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Ourlivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Martellivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Martellivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Picornavirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Picornavirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Reovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Reovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Tolivirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Tolivirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Tymovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Tymovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Nidovirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Nidovirales"
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Nodamuvirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Nodamuvirales"
# ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Stellavirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Stellavirales"

ggt6$data[which(ggt6$data$rdrp_scan_family == regex("Permutotetraviridae", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Permutotetraviridae"

ggt6$data[which(ggt6$data$status == "internal_branch"),]$rdrp_scan_detected_vir_group <- "internal_branch"

ggt6$data$rdrp_scan_assigned <- "yes"
ggt6$data[which(ggt6$data$status == "internal_branch"),]$rdrp_scan_assigned <- "internal_branch"
ggt6$data[which(ggt6$data$status == "novel"),]$rdrp_scan_assigned <- "novel"
ggt6$data[which(is.na(ggt6$data$rdrp_scan_detected_vir_group)),]$rdrp_scan_assigned <- "no"



ggt6 <- ggt6 + geom_tree(mapping = aes(color = rdrp_scan_detected_vir_group, alpha = rdrp_scan_assigned), size = 0.75, layout = "daylight")

ggt6 <- ggt6 + geom_tippoint(aes(size = status, shape = status), alpha = 0.75)


ggt6 <- ggt6 + geom_treescale()


ggt6 <- ggt6 + scale_color_manual(values = c("internal_branch"= "grey80", 

                                             "Pisuviricota"= "#084594",
                                             "Nidovirales" = "#2171b5",
                                             "Durnavirales" = "#4292c6",
                                             "Sobelivirales" = "#6baed6",
                                             "Picornavirales" = "#9ecae1",
                                             
                                             "Lenarviricota" = "#a63603",
                                             "Wolframvirales" = "#e6550d",
                                             "Ourlivirales" = "#fd8d3c",
                                             
                                             "Negarnaviricota" = "#a50f15",
                                             "Monjiviricetes" = "#de2d26",
                                             "Ellioviricetes" = "#fb6a4a",
                                             "Insthoviricetes" = "#fcae91",
                                             
                                             "Duplornaviricota" = "#6a51a3",
                                             "Ghabrivirales" = "#9e9ac8",
                                             "Reovirales" = "#cbc9e2",
                                             
                                             "Kitrinoviricota" = "#005a32",
                                             "Tolivirales" = "#238b45",
                                             "Nodamuvirales" = "#41ab5d",
                                             "Amarillovirales" = "#74c476",
                                             "Tymovirales" = "#a1d99b",
                                             "Martellivirales" = "#c7e9c0",
                                             
                                             "Permutotetraviridae" = "black"
                                             
                                             ))











ggt6 <- ggt6 + scale_alpha_manual(values = c("internal_branch"= 0.25,"yes"= 1,"no" = 0.5, "novel" = 0.5))
ggt6 <- ggt6 + scale_size_manual(values = c("internal_branch"= 0,"reference"= 0,"novel" = 1.5))
ggt6 <- ggt6 + scale_shape_manual(values = c("internal_branch"= NA,"reference"= NA,"novel" = 8))


ggt6 <- ggt6 + theme(legend.position = "bottom", legend.direction = "vertical")
#print(ggt6)

ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_daylight_3.pdf", 
       plot = ggt6, device = "pdf", height = 29, width = 18, units = "cm")



