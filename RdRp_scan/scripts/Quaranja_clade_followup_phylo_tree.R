library(tidyverse)
library(ggtree)
library(treeio)
library(seqinr)
library(lubridate)

################## Tree and metadata 
tree_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/tree/aa_tree/Quaranja_clade_aa_v31_aln_outg.fasta.gap50_clean.fasta.contree"
tree <- read.iqtree(tree_file)
tip_labels <- tree@phylo$tip.label
tip_accessions <- str_extract(tip_labels, '^[[A-Z]]{2}[[\\_]]*?[[:digit:]]+(?=\\_)')

tip_df <- tibble(tip_labels, tip_accessions)

###### Metadata
tip_mdf_brief <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/updated_metadata_with_novel_seq.txt", 
                  delim = "\t")

tip_mdf_brief$novel_seq <- NA

tip_mdf_brief[which(is.na(tip_mdf_brief$nt_accession)),]$novel_seq <- 1

tip_mdf_brief <- tip_mdf_brief[which(tip_mdf_brief$nt_seq_ID %in% tip_df$tip_labels),]

# tip_mdf_brief[which(tip_mdf_brief$novel_seq == 1),] %>% select(nt_seq_ID, Location)


tip_mdf_brief$tip_label_alpha <- "half-transparent"
tip_mdf_brief[which(tip_mdf_brief$novel_seq == 1),]$tip_label_alpha <-  "normal"

tip_mdf_brief <- tip_mdf_brief %>% relocate(nt_seq_ID, .before = nt_accession)





### order genera


tip_mdf_brief$host_genus <- str_extract(tip_mdf_brief$host_taxon, "^[[:alpha:]]+")


sort(unique(tip_mdf_brief$host_genus))
sort(unique(tip_mdf_brief$host_taxon))

tip_mdf_brief$host_genus <- factor(tip_mdf_brief$host_genus, ordered = TRUE,
                                   levels = c('Culex',
                                              'Aedes'))

unique(tip_mdf_brief$host_genus)

unique(tip_mdf_brief$nt_seq_ID)
tip_mdf_brief$var <- tip_mdf_brief$host_genus





sort(unique(tip_mdf_brief$Country))


tip_mdf_brief$Region <- NA

tip_mdf_brief[which(tip_mdf_brief$Country %in% 
                      c("Guadeloupe")),
]$Region <- "Caribbean"

# tip_mdf_brief[which(tip_mdf_brief$Country %in% 
#                       c("Brazil")),
# ]$Region <- "South America"

tip_mdf_brief[which(tip_mdf_brief$Country %in% 
                      c("Sweden", "Greece")),
]$Region <- "Europe"

# tip_mdf_brief[which(tip_mdf_brief$Country %in% 
#                       c("Ghana")),
# ]$Region <- "West Africa"

# tip_mdf_brief[which(tip_mdf_brief$Country %in% 
#                       c("Kenya")),
# ]$Region <- "Southeast Africa"

# tip_mdf_brief[which(tip_mdf_brief$Country %in% 
#                       c()),
# ]$Region <- "Southwest Asia"

tip_mdf_brief[which(tip_mdf_brief$Country %in% 
                      c("China", "Cambodia")),
              ]$Region <- "Southeast Asia"

sort(unique(tip_mdf_brief[which(is.na(tip_mdf_brief$Region)),]$Country))

unique(tip_mdf_brief$Region)


rownames(tip_mdf_brief) <- tip_mdf_brief$nt_seq_ID


genus_matrix <- tip_mdf_brief %>% select(nt_seq_ID, host_genus, var) %>% 
  spread(key = host_genus, value = var, fill = NA) %>% 
  as.data.frame()

rownames(genus_matrix) <- genus_matrix$nt_seq_ID

for (j in 1:ncol(genus_matrix)){
  genus_matrix[, j] <- as.character(genus_matrix[, j])
  genus_matrix[which(is.na(genus_matrix[,j])), j] <- "na"
}

region_df <- tip_mdf_brief %>%
  select(nt_seq_ID, Region)

genus_matrix <-  left_join(genus_matrix, region_df)
rownames(genus_matrix) <- genus_matrix$nt_seq_ID
genus_matrix$nt_seq_ID <- NULL



novel_tip_vec <- tip_mdf_brief[which(is.na(tip_mdf_brief$nt_accession)),]$tip_labels














#### RdRp ORF TREE
tree@phylo$tip.label
tree@treetext

#######################

tree1 <- tree

root <- rootnode(tree1)

edge = tree1@phylo$edge
edge = as.data.frame(edge)
edge = edge[edge[,2] <= Ntip(tree1), ]

df = tree1@data[, c("node", "UFboot")]

### Assign a mock support value to the root/outgroup node, to be able to keep this branch on the tree.
### Replicate sequence branches will also be assigned a mock 100% support

df[which(is.na(df$UFboot)),]$UFboot <- 100

merge(df, edge, by.x="node", by.y="V1", all.x=F, all.y=T) -> df2
df3= df2[, c(3, 2)]
colnames(df3) = colnames(df)
tree1@data = bind_rows(as_tibble(df3), df)

ggt2 <- ggtree(tree1, 
               aes(size=cut(UFboot, c(0, 70, 100))), alpha  = 1/3
) %<+% tip_mdf_brief


ggt2 <- ggt2 + geom_tiplab(data = ggt2$data[which(!is.na(ggt2$data$tip_label_alpha)),], 
  align = TRUE,
  linetype = "dotted",
  linesize = 0.5,
  offset = 0.5, 
  size = 0,
  aes(alpha = tip_label_alpha)
  # size = 2 ### this will put back tip labels
)



####### TEMP
# ggt2 <- ggt2 + geom_label(data = ggt2$data[which(ggt2$data$isTip == F),],
#                           aes(label=node), fill='lightgreen', size = 1.5, alpha = 1/8)
####### 


ggt2 <- ggt2 + geom_treescale(x = 0, y = -5, width = 1, offset = 1, fontsize = 2, linesize = 0.75,color = "grey50")


ggt2 <- ggt2 + scale_alpha_manual(values=c("half-transparent"=1/4,"normal"=1),
                                  guide='legend',
                                  name='tip_label')
ggt2 <- ggt2 + scale_size_manual(values=c('(0,70]'=0.5, '(70,100]'=1),
                                 guide='legend',
                                 name='UFboot')




# virus_clade <- c("Quaranja Charoen-like Quaranjavirus", "Ackless virus",
#                  "Guadeloupe Quaranjavirus", "Adableaps virus", 
#                  "Kimberley Quaranjavirus", "Broome Quaranjavirus 1", "Badu Quaranjavirus",
#                  "Wutai mosquito Quaranjavirus", 
#                  "Fly Quaranjavirus", "Dipteran Quaranjavirus")
# 
# clades_df <- tibble(virus_clade)
# clades_df$node <- NA
# 
# clades_df[which(clades_df$virus_clade == "Quaranja Charoen-like Quaranjavirus"), ]$node <- 86
# clades_df[which(clades_df$virus_clade == "Ackless virus"), ]$node <- 70
# clades_df[which(clades_df$virus_clade == "Guadeloupe Quaranjavirus"), ]$node <- 71
# 
# clades_df[which(clades_df$virus_clade == "Adableaps virus"), ]$node <- 154
# 
# clades_df[which(clades_df$virus_clade == "Kimberley Quaranjavirus"), ]$node <- 76
# clades_df[which(clades_df$virus_clade == "Broome Quaranjavirus 1"), ]$node <- 75
# clades_df[which(clades_df$virus_clade == "Badu Quaranjavirus"), ]$node <- 74
# 
# clades_df[which(clades_df$virus_clade == "Wutai mosquito Quaranjavirus"), ]$node <- 158
# 
# clades_df[which(clades_df$virus_clade == "Fly Quaranjavirus"), ]$node <-80
# clades_df[which(clades_df$virus_clade == "Dipteran Quaranjavirus"), ]$node <- 79
# 
# 
# write.table(x = clades_df,
#             file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/clades.tsv",
#             append = F, sep = "\t", col.names = T,
#             row.names = F, quote = F)

# clades_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/clades.tsv", delim = '\t')
### After first visualization, keep annotation of only the clades that contain new sequences.
### THIS WILL PRODUCE 21 WARNINGS, which is okay
# essential_clades <- c( "AEFV","PCV",
#                        "CxFV","QBV", "CFAV")
# clades_df$virus_clade_essential <- NA
# clades_df[which(clades_df$virus_clade %in% essential_clades), ]$virus_clade_essential <- clades_df[which(clades_df$virus_clade %in% essential_clades), ]$virus_clade
# clades_df$virus_clade_essential <- clades_df$virus_clade
# 
# clades_df$tip <- FALSE
# clades_df[which(clades_df$virus_clade_essential %in% c("Ackless virus",
#                                                        "Guadeloupe Quaranjavirus",
#                                                        "Kimberley Quaranjavirus", "Broome Quaranjavirus 1", "Badu Quaranjavirus",
#                                                        "Fly Quaranjavirus", "Dipteran Quaranjavirus")
#                    ),]$tip <- TRUE
# 
# 
# 
# for (i in 1:nrow(clades_df)){
#   if (clades_df[i, "tip"][[1]]){
#     ggt2 <- ggt2 + geom_cladelab(node=clades_df[i, "node"][[1]],
#                                  label=clades_df[i, "virus_clade_essential"][[1]],
# 
#                                  #align=FALSE, offset = .1, offset.text = 0.1,
#                                  align=FALSE, offset = 0, offset.text = 0,
#                                  textcolour='firebrick3', barcolor='firebrick3',
#                                  fontsize = 2, angle = 0, vjust = 0.4)
# 
# 
#   }else{
#     ggt2 <- ggt2 + geom_cladelab(node=clades_df[i, "node"][[1]],
#                                label=clades_df[i, "virus_clade_essential"][[1]],
# 
#                                #align=FALSE, offset = .1, offset.text = 0.1,
#                                align=FALSE, offset = 0, offset.text = -0.6,
#                                textcolor='firebrick3', barcolor='firebrick3',
#                                barsize=1, fontsize = 2, angle = 0, vjust = 0.4)
# 
#     # ggt2 <- ggt2 + geom_cladelab(node=clades_df[i, "node"][[1]],
#     #                              label=clades_df[i, "virus_clade_essential"][[1]],
#     #
#     #                              #align=FALSE, offset = .1, offset.text = 0.1,
#     #                              align=TRUE, offset = 0, offset.text = -0.45,
#     #                              textcolor='firebrick3', barcolor='firebrick3',
#     #                              barsize=1, fontsize = 0, angle = 10, vjust = 2.5)
#   }
# }

ggt2_hm <- gheatmap(ggt2, genus_matrix, offset=0.55, width=0.12, font.size=2,
                    colnames_angle=-90, hjust=0)



extracted_xlims <- layer_scales(ggt2_hm)$x$range$range

ggt2_hm <- ggt2_hm + geom_text(data = ggt2_hm$data[which(ggt2_hm$data$isTip == T),],
                               aes( x = (extracted_xlims[2] + 0.01), y = y, 
                                    label=Country, 
                                    color = Region), 
                               size=2,
                               hjust = "left")


ggt2_hm <- ggt2_hm + geom_point(data = ggt2_hm$data[which(!is.na(ggt2_hm$data$novel_seq)),],
                          aes( x = (extracted_xlims[2] - 0.32), y = y), 
                          color = "black",
                          size=1, shape = 8)



ggt2_hm <- ggt2_hm + xlim((extracted_xlims[1]), (extracted_xlims[2]+0.1))


#### Coloring of mosquito genus and region
ggt2_hm <- ggt2_hm + scale_fill_manual(breaks=c('Mansonia','Coquillettidia',
                                                'Culex',
                                                'Aedes', 'Haemagogus',
                                                'Diptera',
                                                'Atherigona',
                                                'na',
     
                                                "North America including Central and Caribbean", 
                                                "South America", 
                                                "Europe", 
                                                "West Africa",
                                                "Southeast Africa",
                                                "Southwest Asia",
                                                "Southeast Asia and Australia"                                      
                                                ),
                                       values=c("#44AA99","#117733",
                                                "#332288",
                                                "#882255", "#CC6677",
                                                "grey30",
                                                "grey30",
                                                "grey90",
                                       
                                                "#cbc9e2", 
                                                "#9e9ac8", 
                                                "#0570b0", 
                                                "#c2e699",
                                                "#78c679",
                                                "#238443",
                                                "#d94701"
                                                ))

ggt2_hm <- ggt2_hm + scale_color_manual(breaks=c('Mansonia','Coquillettidia',
                                                 'Culex',
                                                 'Aedes', 'Haemagogus',
                                                 'Diptera',
                                                 'Atherigona',
                                                 'na',
                                                 
                                                 "North America including Central and Caribbean", 
                                                 "South America", 
                                                 "Europe", 
                                                 "West Africa",
                                                 "Southeast Africa",
                                                 "Southwest Asia",
                                                 "Southeast Asia and Australia"                                      
),
values=c("#44AA99","#117733",
         "#332288",
         "#882255", "#CC6677",
         "grey30",
         "grey30",
         "grey90",
         
         "#cbc9e2", 
         "#9e9ac8", 
         "#0570b0", 
         "#c2e699",
         "#78c679",
         "#238443",
         "#d94701"
))





#### add space at the bottom
extracted_ylims <- layer_scales(ggt2)$y$range$range
ggt2_hm <- ggt2_hm + ylim((extracted_ylims[1]-22), (extracted_ylims[2]+5))
ggt2_hm <- ggt2_hm + theme_tree(legend.position="none")
ggt2_hm <- ggt2_hm + theme(plot.margin = unit(c(0, 0, 0, 0), "inch"))

ggt2_hm_legend <- ggt2_hm + theme_tree(legend.position="bottom")
# 
# print(ggt2_hm)
ggt2_hm_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Quaranja_clade_aa_aln.ggtree_hm.pdf"
ggsave(plot = ggt2_hm, filename = ggt2_hm_file, device = "pdf", width = 18, height = 20, units = "cm")

# ggt2_hm_file_legend <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Quaranja_clade_aa_aln.ggtree_hm_legend.pdf"
# ggsave(plot = ggt2_hm_legend, filename = ggt2_hm_file_legend, device = "pdf", width = 50, height = 50, units = "cm")














####################### 
####################### TREES. SUPPLEMENTARY FIGURE
####################### 
tip_mdf_brief$nt_seq_ID
tip_mdf_brief$var <- tip_mdf_brief$host_genus

genus_matrix <- tip_mdf_brief %>% select(nt_seq_ID, host_genus, var) %>% 
  spread(key = host_genus, value = var, fill = NA) %>% 
  as.data.frame()

rownames(genus_matrix) <- genus_matrix$nt_seq_ID
genus_matrix$nt_seq_ID <- NULL

for (j in 1:ncol(genus_matrix)){
  genus_matrix[, j] <- as.character(genus_matrix[, j])
  genus_matrix[which(is.na(genus_matrix[,j])), j] <- "na"
}

novel_tip_vec <- tip_mdf_brief[which(tip_mdf_brief$novel_seq == 1),]$tip_labels

#### nsORF TREE
tree@phylo$tip.label
tree@treetext

tree1 <- tree

root <- rootnode(tree1)

edge = tree1@phylo$edge
edge = as.data.frame(edge)
edge = edge[edge[,2] <= Ntip(tree1), ]

df = tree1@data[, c("node", "UFboot")]

### Assign a mock support value to the root/outgroup node, to be able to keep this branch on the tree.
### Replicate sequence branches will also be assigned a mock 100% support

df[which(is.na(df$UFboot)),]$UFboot <- 100

merge(df, edge, by.x="node", by.y="V1", all.x=F, all.y=T) -> df2
df3= df2[, c(3, 2)]
colnames(df3) = colnames(df)
tree1@data = bind_rows(as_tibble(df3), df)

ggt_sup2 <- ggtree(tree1, 
               aes(size=cut(UFboot, c(0, 70, 100))), alpha  = 1/4
) %<+% tip_mdf_brief
ggt_sup2 <- ggt_sup2 + geom_tiplab(data = ggt2$data[which(!is.na(ggt2$data$tip_label_alpha)),],
  align = TRUE,
  linetype = "dotted",
  linesize = 0.5,
  offset = 0, 
  aes(alpha = tip_label_alpha),
  size = 2 ### this will put back tip labels
  
)

####### TEMP
ggt_sup2 <- ggt_sup2 + geom_label(data = ggt_sup2$data,
                                  # data = ggt_sup2$data[which(ggt_sup2$data$isTip == F),],
                          aes(label=node), fill='lightgreen', size = 1.5, alpha = 1/8)
####### 


ggt_sup2 <- ggt_sup2 + geom_treescale(x = 0, y = -5, width = 1, offset = 1, fontsize = 2, linesize = 0.75)

ggt_sup2 <- ggt_sup2 + scale_alpha_manual(values=c("half-transparent"=1/4,"normal"=1),
                                  guide='legend',
                                  name='tip_label')
ggt_sup2 <- ggt_sup2 + scale_size_manual(values=c('(0,70]'=0.25, '(70,100]'=0.75),
                                 guide='legend',
                                 name='UFboot')

# clades_df <- clades_df_S
# clades_df$node <- NA
# clades_df[which(clades_df$virus_clade == "AEFV"), ]$node <- 270
# clades_df[which(clades_df$virus_clade == "KRV"), ]$node <- 249
# clades_df[which(clades_df$virus_clade == "TrFV"), ]$node <- 143
# clades_df[which(clades_df$virus_clade == "PaRV"), ]$node <- 142
# clades_df[which(clades_df$virus_clade == "AenoFV"), ]$node <- 284
# clades_df[which(clades_df$virus_clade == "HANKV"), ]$node <- 282
# clades_df[which(clades_df$virus_clade == "FalV"), ]$node <- 136
# clades_df[which(clades_df$virus_clade == "OSFV"), ]$node <- 135
# clades_df[which(clades_df$virus_clade == "XFV"), ]$node <- 133
# clades_df[which(clades_df$virus_clade == "MFV"), ]$node <- 134
# clades_df[which(clades_df$virus_clade == "CFAV"), ]$node <- 250
# clades_df[which(clades_df$virus_clade == "MAFV"), ]$node <- 101
# clades_df[which(clades_df$virus_clade == "SbFV"), ]$node <- 100
# clades_df[which(clades_df$virus_clade == "CsFV"), ]$node <- 90
# clades_df[which(clades_df$virus_clade == "MECDV"), ]$node <- 236
# clades_df[which(clades_df$virus_clade == "CLBOV"), ]$node <- 237
# clades_df[which(clades_df$virus_clade == "KRBV"), ]$node <- 231
# clades_df[which(clades_df$virus_clade == "AnFV"), ]$node <- 229
# clades_df[which(clades_df$virus_clade == "McPV"), ]$node <- 83
# clades_df[which(clades_df$virus_clade == "HaCV"), ]$node <- 84
# clades_df[which(clades_df$virus_clade == "NIEV"), ]$node <- 225
# clades_df[which(clades_df$virus_clade == "NAKV"), ]$node <- 76
# clades_df[which(clades_df$virus_clade == "CuCuV"), ]$node <- 77
# clades_df[which(clades_df$virus_clade == "PCV"), ]$node <- 223
# clades_df[which(clades_df$virus_clade == "CTFV"), ]$node <- 208
# clades_df[which(clades_df$virus_clade == "QBV"), ]$node <- 194
# clades_df[which(clades_df$virus_clade == "CxFV"), ]$node <- 151
# 
# write.table(x = clades_df,
#             file = paste0(tree_file, ".clades.tsv"),
#             append = F, sep = "\t", col.names = T,
#             row.names = F, quote = F)

# clades_df <- read_delim(paste0(tree_file, ".clades.tsv"), delim = '\t')
# ### After first visualization, keep annotation of only the clades that contain new sequences.
# ### THIS WILL PRODUCE 21 WARNINGS, which is okay
# essential_clades <- c( "AEFV","PCV",
#                        "CxFV","QBV", "CFAV")
# clades_df$virus_clade_essential <- NA
# clades_df[which(clades_df$virus_clade %in% essential_clades), ]$virus_clade_essential <- clades_df[which(clades_df$virus_clade %in% essential_clades), ]$virus_clade
# 
# 
# 
# 
# 
# 
# for (i in 1:nrow(clades_df)){
#   ggt_sup2 <- ggt_sup2 + geom_cladelab(node=clades_df[i, "node"][[1]],
#                                label=clades_df[i, "virus_clade_essential"][[1]],
#                                #align=FALSE, offset = .1, offset.text = 0.1,
#                                align=TRUE, offset = -0.025, offset.text = -0.4,
#                                textcolor='firebrick3', barcolor='firebrick3',
#                                barsize=0.75, fontsize = 3.5)
# }

ggt_sup2_hm <- gheatmap(ggt_sup2, genus_matrix, offset=1.5, width=0.6, font.size=3,
                    colnames_angle=-90, hjust=0)



#### Coloring of mosquito genus and region
ggt2_hm <- ggt2_hm + scale_fill_manual(breaks=c('Mansonia','Coquillettidia',
                                                'Culex',
                                                'Aedes', 'Haemagogus',
                                                'Diptera',
                                                'Atherigona',
                                                'na',
                                                
                                                "North America including Central and Caribbean", 
                                                "South America", 
                                                "Europe", 
                                                "West Africa",
                                                "Southeast Africa",
                                                "Southwest Asia",
                                                "Southeast Asia and Australia"                                      
),
values=c("#44AA99","#117733",
         "#332288",
         "#882255", "#CC6677",
         "grey30",
         "grey30",
         "grey90",
         
         "#cbc9e2", 
         "#9e9ac8", 
         "#0570b0", 
         "#c2e699",
         "#78c679",
         "#238443",
         "#d94701"
))



#### add space at the bottom
extracted_ylims <- layer_scales(ggt_sup2)$y$range$range
ggt_sup2_hm <- ggt_sup2_hm + ylim((extracted_ylims[1]-22), (extracted_ylims[2]+5))
ggt_sup2_hm <- ggt_sup2_hm + theme_tree(legend.position="none")
#ggt_sup2_hm <- ggt_sup2_hm + theme(plot.margin = unit(c(0, 0, 0, 0), "inch"))


print(ggt_sup2_hm)
ggt_sup2_hm_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Quaranja_clade/Quaranja_clade_aa_aln.ggtree_hm_supplementary.pdf"
ggsave(plot = ggt_sup2_hm, filename = ggt_sup2_hm_file, device = "pdf", width = 21, height = 60, units = "cm")
