library(tidyverse)
library(ggtree)
library(ggrepel)
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
tree <- treeio::read.newick(file = tree_file)
#### Quick tree check
# ggtree(tree)



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



# ggt4_v3 <- ggtree(tree, layout="rectangular")
# save(ggt4_v3,
#      file = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_rectangular_v3.rdata")

load("/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_rectangular_v3.rdata")

tree_tip_df <- tibble(tip_labels_raw)
tree_tip_df1 <- left_join(tree_tip_df, tips_df1, by  = "tip_labels_raw")

tree_tip_df1$phylum <- tree_tip_df1$rdrp_scan_phylum
tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_phylum)),]$phylum <- tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_phylum)),]$BestHit_phylum 

rownames(tree_tip_df1) <- tree_tip_df1$tip_labels_raw





tree_tip_df1$merged_phylum <- tree_tip_df1$rdrp_scan_phylum
tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_phylum)),]$merged_phylum <- tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_phylum)),]$BestHit_phylum

tree_tip_df1$merged_class <- tree_tip_df1$rdrp_scan_class
tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_class)),]$merged_class <- tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_class)),]$BestHit_class

tree_tip_df1$merged_order <- tree_tip_df1$rdrp_scan_order
tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_order)),]$merged_order <- tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_order)),]$BestHit_order

tree_tip_df1$merged_family <- tree_tip_df1$rdrp_scan_family
tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_family)),]$merged_family <- tree_tip_df1[which(!is.na(tree_tip_df1$BestHit_family)),]$BestHit_family

tree_tip_df1$merged_genus<- tree_tip_df1$BestHit_genus

tree_tip_df1$merged_vir_sp <- str_remove_all(
  str_remove_all(
    str_extract(tree_tip_df1$tip_labels, "\\_\\_[[:alnum:]\\_\\.]+$"), 
                                                            "^\\_\\_"), 
                                             "\\_$" )
tree_tip_df1[which(!is.na(tree_tip_df1$suggested_vir_sp)),]$merged_vir_sp <- tree_tip_df1[which(!is.na(tree_tip_df1$suggested_vir_sp)),]$suggested_vir_sp
tree_tip_df1[which(is.na(tree_tip_df1$merged_vir_sp) & 
                     str_detect(tree_tip_df1$tip_labels, "virus")),]$merged_vir_sp <- tree_tip_df1[which(is.na(tree_tip_df1$merged_vir_sp) & 
                                                                                                                                             str_detect(tree_tip_df1$tip_labels, "virus")),]$tip_labels

tree_tip_df1$merged_vir_sp <-   str_remove_all(str_remove_all(str_remove_all(tree_tip_df1$merged_vir_sp, "[[:print:]]+\\|"), 
                                                              "^[\\_]+"), 
                                               "[\\_]+$" )


for(i in 1:nrow(tree_tip_df1)){
  taxa_annotated <- unique(as.vector(na.omit(c(tree_tip_df1[i,]$rdrp_scan_phylum, 
                              tree_tip_df1[i,]$rdrp_scan_class, 
                              tree_tip_df1[i,]$rdrp_scan_order, 
                              tree_tip_df1[i,]$rdrp_scan_family))))
  
  if (length(taxa_annotated)){
    for (t in taxa_annotated){
      if(is.na(str_detect(tree_tip_df1[i,]$merged_vir_sp, t)) == FALSE ) {
        # print(t)
        # print(tree_tip_df1[i,]$merged_vir_sp)
        tree_tip_df1[i,]$merged_vir_sp <- str_remove_all(str_remove_all(tree_tip_df1[i,]$merged_vir_sp, t), 
                                                         "^[\\_]+")
        # print(tree_tip_df1[i,]$merged_vir_sp)
      }
    }
  }
}



tree_tip_df1$merged_vir_sp <- str_remove_all(str_remove_all(tree_tip_df1$merged_vir_sp, "NOT_ASSIGNED"), 
                                                            "^[\\_]+")







### check order
sum(rownames(tree_tip_df1) == tip_labels_raw) == nrow(tree_tip_df1) 






############# PLOT
######## Figure for manual clade annotation
################################################################################
################################################################################
################################################################################

### Import taxonomy based on DIAMOND blastp and RdRp-scan. Taxonomic gaps for higher ranks were filled out in ./viz_prevalence.R
df_tax <- read_delim("/full_path_to/wd/RdRp_scan/analysis/read_mapping/best_hit_taxonomy.txt", delim = "\t")
df_tax1 <- df_tax %>% select(orf_id, pident, evalue, sseqid_BE = sseqid, rdrp_scan_evalue,
                             phylum_BE = phylum, class_BE = class, order_BE = order, 
                             family_BE = family, genus_BE = genus, species_BE = species, 
                             suggested_vir_sp)

df_tax1$label <- paste0("'", df_tax1$orf_id, "'") #### sometimes raw tip labels on the tree are read with single quotes, not sure why, seems to be different on different computers, perhaps versions of packages
### this seems to happen only if there are special characters in the tip label besides underscore
### All novel sequences have \\. in the name so, I can just add '' to all orf_id names for now

# df_tax1$label <- df_tax1$orf_id  ### alternative





tree_tip_df2 <- tree_tip_df1 %>% select(tip_labels_raw)
df_tax2 <- left_join(tree_tip_df1, df_tax1, by = c("tip_labels_raw"="label")) %>%
  unclass() %>%
  as.data.frame(stringsAsFactors=TRUE)
rownames(df_tax2) <- df_tax2$tip_labels_raw
hm_matrix <- df_tax2 %>% select(pident) %>% 
  as.data.frame()


################################################################################

ggt41 <- ggt4_v3  
ggt41$layers[[1]]$geom$default_aes$size <- 0
ggt41$layers[[1]]$geom$default_aes$alpha <- 0

ggt6 <- ggt41 %<+% tree_tip_df1

save(ggt6,
     file = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_rectangular_v3_data.rdata")




##### collapse nodes
#### select nodes to collapse
nodes_to_collapse <- sort(unique(c(4134, 3799, 3874, 3915, 3945,
                       4036, 4058, 4096, 4110,
                       4122, 4535, 4518, 4526, 4561,
                       4568, 4580, 4083, 4582,
                       4646, 4655, 4714, 4723, 
                       4790, 4811, 4847, 4979,
                       5045, 5050, 5083, 5117,
                       5135, 5156, 5163, 5177,
                       5188, 5194, 5197, 5199, 
                       5221, 5223, 5257, 5251,
                       5372, 5492, 5587, 5592,
                       5611, 5625, 5639, 5655,
                       5659, 5662, 5675, 5679,
                       5687, 5691, 5750, 5798,
                       5805, 5815, 5902, 6070, 6290,
                       6416, 6426, 6629, 7035,
                       7040, 7087, 7108, 7113,
                       7227, 7252, 4405)))




for (nd in nodes_to_collapse){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )

}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse)), 
                           shape=23, size=2, 
                           fill = "grey80")




nodes_to_collapse_Tolivirales <- c(3713, 3731,3816, 4091, 4771)
for (nd in nodes_to_collapse_Tolivirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Tolivirales)), 
                           shape=23, size=2, 
                           fill = "#238b45")
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Tolivirales),
                              label = "Tolivirales"),
                          size=2, 
                          color = "#238b45", nudge_x = 0.5
)



nodes_to_collapse_Amarillovirales <- c(4170, 4816)
for (nd in nodes_to_collapse_Amarillovirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Amarillovirales)), 
                           shape=23, size=2, 
                           fill = "#74c476")
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Amarillovirales),
                              label = "Amarillovirales"),
                          size=2, 
                          color = "#74c476", nudge_x = 0.5
)

nodes_to_collapse_Martellivirales <- c(4231, 4246, 4369)
for (nd in nodes_to_collapse_Martellivirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Martellivirales)), 
                           shape=23, size=2, 
                           fill = "#c7e9c0")
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Martellivirales),
                              label = "Martellivirales"),
                          size=2, 
                          color = "#c7e9c0", nudge_x = 0.5
)




nodes_to_collapse_Hepelivirales <- c(4310)
for (nd in nodes_to_collapse_Hepelivirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Hepelivirales)), 
                           shape=23, size=2, 
                           fill = "#005a32" #### color for Kitrinoviricota because I don't annotate at the order level Hapelivirales if there's no novel viruses there to save on color shades used
                           )
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Hepelivirales),
                              label = "Hepelivirales"),
                          size=2, 
                          color = "#005a32", #### color for Kitrinoviricota because I don't annotate at the order level Hapelivirales if there's no novel viruses there to save on color shades used
                          nudge_x = 0.5
)




nodes_to_collapse_Tymovirales <- c(4396, 4388)
for (nd in nodes_to_collapse_Tymovirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Tymovirales)), 
                           shape=23, size=2, 
                           fill = "#a1d99b" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Tymovirales),
                              label = "Tymovirales"),
                          size=2, 
                          color = "#a1d99b", nudge_x = 0.5
)

nodes_to_collapse_Nodamuvirales <- c(4419)
for (nd in nodes_to_collapse_Nodamuvirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Nodamuvirales)), 
                           shape=23, size=2, 
                           fill = "#41ab5d" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Nodamuvirales),
                              label = "Nodamuvirales"),
                         size=2, 
                         color = "#41ab5d", nudge_x = 0.5
)




nodes_to_collapse_Picornavirales <- c(4511, 4867, 5130, 5125, 5226, 
                                      5242, 5183, 5205, 5217, 5280,
                                      5359, 5419, 5433, 5445, 5464,
                                      5495, 5500, 5540, 5581, 5599
                                      )
for (nd in nodes_to_collapse_Picornavirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Picornavirales)), 
                           shape=23, size=2, 
                           fill = "#a6bddb" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Picornavirales),
                              label = "Picornavirales"),
                          size=2, 
                          color = "#a6bddb", nudge_x = 0.5
)




nodes_to_collapse_Patatavirales <- c(5666)
for (nd in nodes_to_collapse_Patatavirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Patatavirales)), 
                           shape=23, size=2, 
                           fill = "#034e7b" ## Pisuviricota color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Patatavirales),
                              label = "Patatavirales"),
                          size=2, 
                          color = "#034e7b", ## Pisuviricota color
                          nudge_x = 0.5
)




nodes_to_collapse_Sobelivirales <- c(4609, 4760)
for (nd in nodes_to_collapse_Sobelivirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Sobelivirales)), 
                           shape=23, size=2, 
                           fill = "#74a9cf" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Sobelivirales),
                              label = "Sobelivirales"),
                          size=2, 
                          color = "#74a9cf", nudge_x = 0.5
)




nodes_to_collapse_Hypovirus <- c(5703)
for (nd in nodes_to_collapse_Hypovirus){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Hypovirus)), 
                           shape=23, size=2, 
                           fill = "#3690c0" ## Durnavirales color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Hypovirus),
                              label = "Hypovirus"),
                          size=2, 
                          color = "#3690c0", ## Durnavirales color
                          nudge_x = 0.5
)



nodes_to_collapse_Astrovirus <- c(5721)
for (nd in nodes_to_collapse_Astrovirus){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Astrovirus)), 
                           shape=23, size=2, 
                           fill = "#d0d1e6" ## Stellavirales color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Astrovirus),
                              label = "Astrovirus"),
                          size=2, 
                          color = "#d0d1e6", ## Stellavirales color
                          nudge_x = 0.5
)



nodes_to_collapse_Narna_like <- c(5781, 6050)
for (nd in nodes_to_collapse_Narna_like){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Narna_like)), 
                           shape=23, size=2, 
                           fill = "#e6550d" ## Wolframvirales color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Narna_like),
                              label = "Narna_like"),
                          size=2, 
                          color = "#e6550d", ## Wolframvirales color
                          nudge_x = 0.5
)


nodes_to_collapse_Wolfram_Ourli_virales <- c(5858)
for (nd in nodes_to_collapse_Wolfram_Ourli_virales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Wolfram_Ourli_virales)), 
                           shape=23, size=2, 
                           fill = "#a63603" ## Lenaviricota color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Wolfram_Ourli_virales),
                              label = "Wolfram_Ourli_virales"),
                          size=2, 
                          color = "#a63603", ## Lenaviricota color
                          nudge_x = 0.5
)



nodes_to_collapse_Cryppa_Wolfram_virales <- c(6103)
for (nd in nodes_to_collapse_Cryppa_Wolfram_virales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Cryppa_Wolfram_virales)), 
                           shape=23, size=2, 
                           fill = "#a63603" ## Lenaviricota color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Cryppa_Wolfram_virales),
                              label = "Cryppa_Wolfram_virales"),
                          size=2, 
                          color = "#a63603", ## Lenaviricota color
                          nudge_x = 0.5
)


nodes_to_collapse_Levivirales <- c(6298)
for (nd in nodes_to_collapse_Levivirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Levivirales)), 
                           shape=23, size=2, 
                           fill = "#a63603" ## Lenaviricota color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Levivirales),
                              label = "Levivirales"),
                          size=2, 
                          color = "#a63603", ## Lenaviricota color
                          nudge_x = 0.5
)


nodes_to_collapse_Chunqiuviricetes <- c(6764)
for (nd in nodes_to_collapse_Chunqiuviricetes){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Chunqiuviricetes)), 
                           shape=23, size=2, 
                           fill = "#a50f15" ## Negarnaviricota color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Chunqiuviricetes),
                              label = "Chunqiuviricetes"),
                          size=2, 
                          color = "#a50f15", ## Negarnaviricotacolor
                          nudge_x = 0.5
)

nodes_to_collapse_Milneviricetes <- c(6755)
for (nd in nodes_to_collapse_Milneviricetes){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Milneviricetes)), 
                           shape=23, size=2, 
                           fill = "#a50f15" ## Negarnaviricota color
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Milneviricetes),
                              label = "Milneviricetes"),
                          size=2, 
                          color = "#a50f15", ## Negarnaviricotacolor
                          nudge_x = 0.5
)

nodes_to_collapse_Monjiviricetes <- c(6732)
for (nd in nodes_to_collapse_Monjiviricetes){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Monjiviricetes)), 
                           shape=23, size=2, 
                           fill = "#de2d26" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Monjiviricetes),
                              label = "Monjiviricetes"),
                          size=2, 
                          color = "#de2d26", 
                          nudge_x = 0.5
)


nodes_to_collapse_Ellioviricetes <- c(6623)
for (nd in nodes_to_collapse_Ellioviricetes){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Ellioviricetes)), 
                           shape=23, size=2, 
                           fill = "#fb6a4a" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Ellioviricetes),
                              label = "Ellioviricetes"),
                          size=2, 
                          color = "#fb6a4a", 
                          nudge_x = 0.5
)


nodes_to_collapse_Ghabrivirales <- c(6865, 6946, 6920, 6818)
for (nd in nodes_to_collapse_Ghabrivirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Ghabrivirales)), 
                           shape=23, size=2, 
                           fill = "#8856a7" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Ghabrivirales),
                              label = "Ghabrivirales"),
                          size=2, 
                          color = "#8856a7", 
                          nudge_x = 0.5
)




nodes_to_collapse_Birnaviridae <- c(7099)
for (nd in nodes_to_collapse_Birnaviridae){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Birnaviridae)), 
                           shape=23, size=2, 
                           fill = "pink" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Birnaviridae),
                              label = "Birnaviridae"),
                          size=2, 
                          color = "pink", 
                          nudge_x = 0.5
)




nodes_to_collapse_Mindivirales <- c(7027)
for (nd in nodes_to_collapse_Mindivirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Mindivirales)), 
                           shape=23, size=2, 
                           fill = "#810f7c" ## Duplornaviricota
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Mindivirales),
                              label = "Mindivirales"),
                          size=2, 
                          color = "#810f7c", ## Duplornaviricota
                          nudge_x = 0.5
)




nodes_to_collapse_Reovirales <- c(6977)
for (nd in nodes_to_collapse_Reovirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Reovirales)), 
                           shape=23, size=2, 
                           fill = "#8c96c6" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Reovirales),
                              label = "Reovirales"),
                          size=2, 
                          color = "#8c96c6", 
                          nudge_x = 0.5
)

nodes_to_collapse_Durnavirales <- c(7163,7168,7178,7218,7258)
for (nd in nodes_to_collapse_Durnavirales){
  tryCatch(
    expr = {
      ggt6 <- ggt6 %>% collapse(node=nd)
    },
    warning = function(w){
      print(nd)
    }
  )
  
}
ggt6 <- ggt6 + geom_point2(aes(subset=(node %in% nodes_to_collapse_Durnavirales)), 
                           shape=23, size=2, 
                           fill = "#3690c0" 
)
ggt6 <- ggt6 + geom_text2(aes(subset=(node %in% nodes_to_collapse_Durnavirales),
                              label = "Durnavirales"),
                          size=2, 
                          color = "#3690c0", 
                          nudge_x = 0.5
)



























ggt6$data[which(is.na(ggt6$data$status)),]$status <- "internal_branch"
ggt6$data[which(is.na(ggt6$data$status)),]$status <- factor(ggt6$data[which(is.na(ggt6$data$status)),]$status, 
                                                            levels = c("internal_branch", "reference", "novel" ), ordered = T)


####### Temporary node label for annotation
# ggt6 <- ggt6 + geom_label(data = ggt6$data[which(ggt6$data$isTip == F),], 
#                           aes(label=node), fill='lightgreen', size = 1.5, alpha = 1/8)













#### only taxa for which novel sequences were detected
ggt6$data$rdrp_scan_detected_vir_group <- NA
#### phylum
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Pisuviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Pisuviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Kitrinoviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Kitrinoviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Duplornaviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Duplornaviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Lenarviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Lenarviricota"
ggt6$data[which(ggt6$data$rdrp_scan_phylum == regex("Negarnaviricota", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Negarnaviricota"

####
# temp_df <- ggt6$data[which(!ggt6$data$rdrp_scan_phylum %in% c("Pisuviricota", "Kitrinoviricota", "Duplornaviricota", "Lenarviricota", "Negarnaviricota") & !ggt6$data$status == "internal_branch"),]


#### class
ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Monjiviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Monjiviricetes"
ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Ellioviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Ellioviricetes"
ggt6$data[which(ggt6$data$rdrp_scan_class == regex("Insthoviricetes", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Insthoviricetes"
#### order
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
ggt6$data[which(ggt6$data$rdrp_scan_order == regex("Stellavirales", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Stellavirales"

#### family
ggt6$data[which(ggt6$data$rdrp_scan_family == regex("Permutotetraviridae", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Permutotetraviridae"
ggt6$data[which(ggt6$data$rdrp_scan_family == regex("Birnaviridae", ignore_case = T)),]$rdrp_scan_detected_vir_group <- "Birnaviridae"



ggt6$data[which(ggt6$data$status == "internal_branch"),]$rdrp_scan_detected_vir_group <- "internal_branch"
ggt6$data$rdrp_scan_assigned <- "yes"
ggt6$data[which(is.na(ggt6$data$rdrp_scan_detected_vir_group)),]$rdrp_scan_assigned <- "no"
ggt6$data[which(ggt6$data$status == "internal_branch"),]$rdrp_scan_assigned <- "internal_branch"
ggt6$data[which(ggt6$data$status == "novel"),]$rdrp_scan_assigned <- "novel"


ggt6 <- ggt6 + geom_tree(mapping = aes(color = rdrp_scan_detected_vir_group, alpha = rdrp_scan_assigned, 
                                       linetype = rdrp_scan_assigned),
                         size = 0.5, layout = "rectangular")

ggt6 <- ggt6 + geom_tippoint(aes(size = status, shape = status), alpha = 0.75)


####### add tip annotation, just keep the alignment line
ggt6 <- ggt6 + geom_tiplab(
  align = TRUE,
  linetype = "dashed",
  linesize = 0.1,
  offset = 0,
  size = 0, ### just keep the line
  aes(color = rdrp_scan_detected_vir_group, 
      alpha = rdrp_scan_assigned)
)

####### Additional annotation
ggt6 <- ggt6 + geom_text(data = ggt6$data[which(ggt6$data$isTip == T),],
                         aes( x = (8.023116 - 0.7), y = y, 
                              label=merged_vir_sp, 
                              color = rdrp_scan_detected_vir_group, 
                              alpha = rdrp_scan_assigned), 
                         size=0.5,
                         hjust = "left")
ggt6 <- ggt6 + geom_text(data = ggt6$data[which(ggt6$data$isTip == T),],
                         aes( x = (8.023116 - 0.4), y = y, 
                              label=merged_family, 
                              color = rdrp_scan_detected_vir_group, 
                              alpha = rdrp_scan_assigned), 
                         size=0.5,
                         hjust = "left")
ggt6 <- ggt6 + geom_text(data = ggt6$data[which(ggt6$data$isTip == T),],
                         aes( x = (8.023116 - 0.2), y = y, 
                              label=merged_order, 
                              color = rdrp_scan_detected_vir_group, 
                              alpha = rdrp_scan_assigned), 
                         size=0.5,
                         hjust = "left")
ggt6 <- ggt6 + geom_text(data = ggt6$data[which(ggt6$data$isTip == T),],
                         aes( x = (8.023116 + 0), y = y, 
                              label=merged_class, 
                              color = rdrp_scan_detected_vir_group, 
                              alpha = rdrp_scan_assigned), 
                         size=0.5,
                         hjust = "left")
ggt6 <- ggt6 + geom_text(data = ggt6$data[which(ggt6$data$isTip == T),],
                         aes( x = (8.023116 + 0.2), y = y, 
                              label=merged_phylum, 
                              color = rdrp_scan_detected_vir_group, 
                              alpha = rdrp_scan_assigned), 
                         size=0.5,
                         hjust = "left")
# ggt6 <- ggt6 + geom_text(data = ggt6$data[which(ggt6$data$isTip == T),],
#                          aes( x = (8.023116 + 1.5), y = y, 
#                               label=tip_labels, 
#                               color = rdrp_scan_detected_vir_group, 
#                               alpha = rdrp_scan_assigned), 
#                          size = 0.5,
#                          hjust = "left")







#### add pident heatmap
ggt6  <- gheatmap(ggt6, hm_matrix, offset= -0.1 , width=0.05, font.size=0,
                  colnames_angle=0, hjust=0)




ggt6 <- ggt6 + geom_treescale(x = 0, y = 1450)


ggt6 <- ggt6 + scale_color_manual(values = c("internal_branch"= "grey80", 
                                             
                                             "Pisuviricota"= "#034e7b",
                                             "Nidovirales" = "#0570b0",
                                             "Durnavirales" = "#3690c0",
                                             "Sobelivirales" = "#74a9cf",
                                             "Picornavirales" = "#a6bddb",
                                             "Stellavirales" = "#d0d1e6",
                                             
                                             "Lenarviricota" = "#a63603",
                                             "Wolframvirales" = "#e6550d",
                                             "Ourlivirales" = "#fd8d3c",
                                             
                                             "Negarnaviricota" = "#a50f15",
                                             "Monjiviricetes" = "#de2d26",
                                             "Ellioviricetes" = "#fb6a4a",
                                             "Insthoviricetes" = "#fcae91",
                                             
                                             "Duplornaviricota" = "#810f7c",
                                             "Ghabrivirales" = "#8856a7",
                                             "Reovirales" = "#8c96c6",
                                             
                                             "Kitrinoviricota" = "#005a32",
                                             "Tolivirales" = "#238b45",
                                             "Nodamuvirales" = "#41ab5d",
                                             "Amarillovirales" = "#74c476",
                                             "Tymovirales" = "#a1d99b",
                                             "Martellivirales" = "#c7e9c0",
                                             
                                             "Permutotetraviridae" = "black",
                                             "Birnaviridae" = "pink"
                                             
))

ggt6 <- ggt6 + scale_alpha_manual(values = c("internal_branch"= 0.25,"yes"= 1,"no" = 0.5, "novel" = 0.5))
ggt6 <- ggt6 + scale_size_manual(values = c("internal_branch"= 0,"reference"= 0,"novel" = 1.5))
ggt6 <- ggt6 + scale_shape_manual(values = c("internal_branch"= NA,"reference"= NA,"novel" = 8))
ggt6 <- ggt6 + scale_linetype_manual(values = c("internal_branch"= "solid", "yes"= "solid", "no" = "dotted", "novel" = "solid"))


#### add space at the bottom
extracted_xlims <- layer_scales(ggt6)$x$range$range
ggt6 <- ggt6 + xlim((extracted_xlims[1]), (extracted_xlims[2]+0.5))


#### legend
# ggt6 <- ggt6 + theme(legend.position = "bottom", legend.direction = "vertical")
ggt6 <- ggt6 + theme(legend.position = "none")


ggsave(filename = "/full_path_to/wd/RdRp_scan/analysis/phylo/R_232orfs_aa_fasttree_tips_clean.ggtree_midpoint_by_figtree_rectangular_4.pdf", 
       plot = ggt6, device = "pdf", height = 42, width = 29.7, units = "cm", limitsize = F)
















