library(tidyverse)
library(ggtree)
library(treeio)
library(seqinr)
library(lubridate)

################## Tree and metadata 
tree_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/tree/aa_tree/Sobeli_clade_aa_v2311_dedup_aln_outg.fasta.gap50_clean.fasta.contree.manual_reroot.txt"
tree <- read.iqtree(tree_file)
tip_labels <- tree@phylo$tip.label
tip_accessions <- str_extract(tip_labels, '^[[A-Z]]{2}[[\\_]]*?[[:digit:]]+(?=\\_)')

tip_df <- tibble(tip_labels, tip_accessions)

###### Metadata
tip_mdf_brief <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/updated_metadata_with_novel_seq.txt", 
                  delim = "\t")

tip_mdf_brief$novel_seq <- NA

tip_mdf_brief[which(is.na(tip_mdf_brief$nt_accession)),]$novel_seq <- 1

tip_mdf_brief <- tip_mdf_brief[which(tip_mdf_brief$nt_seq_ID %in% tip_df$tip_labels),]

# tip_mdf_brief[which(tip_mdf_brief$novel_seq == 1),] %>% select(nt_seq_ID, Location)


tip_mdf_brief$tip_label_alpha <- "half-transparent"
tip_mdf_brief[which(tip_mdf_brief$novel_seq == 1),]$tip_label_alpha <-  "normal"

tip_mdf_brief <- tip_mdf_brief %>% relocate(nt_seq_ID, .before = nt_accession)




#### fill in empty host metadata
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MZ556263"),]$host_taxon <- "Culicidae"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MW434863"),]$host_taxon <- "Culex tarsalis"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MW434872"),]$host_taxon <- "Culex tarsalis"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MW434875"),]$host_taxon <- "Culex tarsalis"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MW434902"),]$host_taxon <- "Aedes aegypti"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MT896210"),]$host_taxon <- "Aedes aegypti"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MT896211"),]$host_taxon <- "Aedes aegypti"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MT361053"),]$host_taxon <- "Aedes albopictus"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MT361055"),]$host_taxon <- "Aedes albopictus"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MT361059"),]$host_taxon <- "Aedes albopictus"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MT361061"),]$host_taxon <- "Aedes albopictus"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MT096519"),]$host_taxon <- "Arthropoda"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MN661097"),]$host_taxon <- "Coquillettidia venezuelensis"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MN661099"),]$host_taxon <- "Coquillettidia venezuelensis"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MN661107"),]$host_taxon <- "Coquillettidia venezuelensis"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MK440646"),]$host_taxon <- "Culex"
tip_mdf_brief[which(tip_mdf_brief$nt_accession == "MH703049"),]$host_taxon <- "Aedes vexans nipponii"


tip_mdf_brief[which(is.na(tip_mdf_brief$host_taxon)),]$nt_accession ## check


### order genera
tip_mdf_brief$host_genus <- str_extract(tip_mdf_brief$host_taxon, "^[[:alpha:]]+")


sort(unique(tip_mdf_brief$host_genus))
sort(unique(tip_mdf_brief$host_taxon))


#### Assign higher level host taxa
tip_mdf_brief$host_taxon1 <- tip_mdf_brief$host_genus

sort(unique(tip_mdf_brief$host_taxon1))










tip_mdf_brief$host_taxon1 <- factor(tip_mdf_brief$host_taxon1, ordered = TRUE,
                                   levels = c('Anopheles', 
                                              'Armigeres',
                                              'Coquillettidia',
                                              'Culiseta',
                                              'Culex',
                                              'Aedes',
                                              'Culicidae',
                                              'Arthropoda'))

unique(tip_mdf_brief$host_taxon1)

unique(tip_mdf_brief$nt_seq_ID)
tip_mdf_brief$var <- tip_mdf_brief$host_taxon1





sort(unique(tip_mdf_brief$Country))

tip_mdf_brief$Region <- NA

tip_mdf_brief[which(tip_mdf_brief$Country %in% 
                      c("USA", "Mexico", "Guadeloupe")),
]$Region <- "North America and Caribbean"

tip_mdf_brief[which(tip_mdf_brief$Country %in%
                      c("Colombia", "Brazil")),
]$Region <- "South America"

tip_mdf_brief[which(tip_mdf_brief$Country %in% 
                      c("Sweden","Spain")),
]$Region <- "Europe"

# tip_mdf_brief[which(tip_mdf_brief$Country %in% 
#                       c("Ghana")),
# ]$Region <- "West Africa"

tip_mdf_brief[which(tip_mdf_brief$Country %in%
                      c("Kenya")),
]$Region <- "Southeast Africa"

# tip_mdf_brief[which(tip_mdf_brief$Country %in% 
#                       c()),
# ]$Region <- "Southwest Asia"

tip_mdf_brief[which(tip_mdf_brief$Country %in% 
                      c("China", "Cambodia","South Korea", "Japan")),
              ]$Region <- "Southeast Asia"

sort(unique(tip_mdf_brief[which(is.na(tip_mdf_brief$Region)),]$Country))


unique(tip_mdf_brief$Region)


rownames(tip_mdf_brief) <- tip_mdf_brief$nt_seq_ID


genus_matrix <- tip_mdf_brief %>% select(nt_seq_ID, host_taxon1, var) %>% 
  spread(key = host_taxon1, value = var, fill = NA) %>% 
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
  offset = 0.2, 
  size = 0,
  aes(alpha = tip_label_alpha)
  # size = 2 ### this will put back tip labels
)



####### TEMP
# ggt2 <- ggt2 + geom_label(data = ggt2$data[which(ggt2$data$isTip == F),],
#                           aes(label=node), fill='lightgreen', size = 1.5, alpha = 1/8)
####### 


ggt2 <- ggt2 + geom_treescale(x = 0, y = -5, width = 0.2, offset = 1, fontsize = 2, linesize = 0.75, color = "grey50")


ggt2 <- ggt2 + scale_alpha_manual(values=c("half-transparent"=1/4,"normal"=1),
                                  guide='legend',
                                  name='tip_label')
ggt2 <- ggt2 + scale_size_manual(values=c('(0,70]'=0.5, '(70,100]'=1),
                                 guide='legend',
                                 name='UFboot')




virus_clade <- c("Comfornsto virus", "Kellev virus",
                 "Kisumu mosquito virus",
                 "Hubei mosquito virus 2 NC_033305",
                 "Marma virus", "Hubei mosquito virus 2",
                 "Culex inatomii luteo-like virus",
                 "Peninsere virus", "Pirounced virus", 
                 "Guadeloupe mosquito virus",
                 "Kwale mosquito virus", "Renna virus",
                 "Colokingly virus", "Aedes sobemo-like virus",
                 "Atrato Sobemo-like virus 5", "Guangzhou sobemo-like virus",
                 "Sichuan mosquito sobemo-like virus",
                 "Wenzhou sobemo-like virus 4 MT096519", 
                 "Wenzhou sobemo-like virus 4 NC_033138")

clades_df <- tibble(virus_clade)
clades_df$node <- NA

clades_df[which(clades_df$virus_clade == "Comfornsto virus"), ]$node <- 143
clades_df[which(clades_df$virus_clade == "Kellev virus"), ]$node <- 81
clades_df[which(clades_df$virus_clade == "Kisumu mosquito virus"), ]$node <- 34
clades_df[which(clades_df$virus_clade == "Hubei mosquito virus 2 NC_033305"), ]$node <- 33
clades_df[which(clades_df$virus_clade == "Marma virus"), ]$node <- 114
clades_df[which(clades_df$virus_clade == "Hubei mosquito virus 2"), ]$node <- 88
clades_df[which(clades_df$virus_clade == "Culex inatomii luteo-like virus"), ]$node <- 13
clades_df[which(clades_df$virus_clade == "Peninsere virus"), ]$node <- 111
clades_df[which(clades_df$virus_clade == "Pirounced virus"), ]$node <- 92
clades_df[which(clades_df$virus_clade == "Guadeloupe mosquito virus"), ]$node <- 145
clades_df[which(clades_df$virus_clade == "Kwale mosquito virus"), ]$node <- 64
clades_df[which(clades_df$virus_clade == "Renna virus"), ]$node <- 62
clades_df[which(clades_df$virus_clade == "Colokingly virus"), ]$node <- 153
clades_df[which(clades_df$virus_clade == "Aedes sobemo-like virus"), ]$node <- 124
clades_df[which(clades_df$virus_clade == "Atrato Sobemo-like virus 5"), ]$node <- 122
clades_df[which(clades_df$virus_clade == "Guangzhou sobemo-like virus"), ]$node <- 125
clades_df[which(clades_df$virus_clade == "Sichuan mosquito sobemo-like virus"), ]$node <- 57
clades_df[which(clades_df$virus_clade == "Wenzhou sobemo-like virus 4 MT096519"), ]$node <- 56
clades_df[which(clades_df$virus_clade == "Wenzhou sobemo-like virus 4 NC_033138"), ]$node <- 58


write.table(x = clades_df,
            file = "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/clades.tsv",
            append = F, sep = "\t", col.names = T,
            row.names = F, quote = F)

clades_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/clades.tsv", delim = '\t')

clades_df$virus_clade_essential <- clades_df$virus_clade

clades_df$tip <- FALSE
clades_df[which(clades_df$virus_clade_essential %in% c("Kellev virus",
                                                       "Kisumu mosquito virus",
                                                       "Hubei mosquito virus 2 NC_033305",
                                                       "Culex inatomii luteo-like virus", 
                                                       "Kwale mosquito virus",
                                                       "Renna virus",
                                                       "Sichuan mosquito sobemo-like virus",
                                                       "Wenzhou sobemo-like virus 4 MT096519", 
                                                       "Wenzhou sobemo-like virus 4 NC_033138")),]$tip <- TRUE



for (i in 1:nrow(clades_df)){
  if (clades_df[i, "tip"][[1]]){
    ggt2 <- ggt2 + geom_cladelab(node=clades_df[i, "node"][[1]],
                                 label=clades_df[i, "virus_clade_essential"][[1]],

                                 #align=FALSE, offset = .1, offset.text = 0.1,
                                 align=FALSE, offset = 0, offset.text = 0,
                                 textcolour='firebrick3', barcolor='firebrick3',
                                 fontsize = 2, angle = 0, vjust = 0.4)


  }else{
    ggt2 <- ggt2 + geom_cladelab(node=clades_df[i, "node"][[1]],
                               label=clades_df[i, "virus_clade_essential"][[1]],

                               #align=FALSE, offset = .1, offset.text = 0.1,
                               align=FALSE, offset = 0, offset.text = -0.1,
                               textcolor='firebrick3', barcolor='firebrick3',
                               barsize=1, fontsize = 2, angle = 0, vjust = 0.4)

    # ggt2 <- ggt2 + geom_cladelab(node=clades_df[i, "node"][[1]],
    #                              label=clades_df[i, "virus_clade_essential"][[1]],
    #
    #                              #align=FALSE, offset = .1, offset.text = 0.1,
    #                              align=TRUE, offset = 0, offset.text = -0.45,
    #                              textcolor='firebrick3', barcolor='firebrick3',
    #                              barsize=1, fontsize = 0, angle = 10, vjust = 2.5)
  }
}

ggt2_hm <- gheatmap(ggt2, genus_matrix, offset=0.225, width=0.325, font.size=2,
                    colnames_angle=-90, hjust=0)



extracted_xlims <- layer_scales(ggt2_hm)$x$range$range

ggt2_hm <- ggt2_hm + geom_text(data = ggt2_hm$data[which(ggt2_hm$data$isTip == T),],
                               aes( x = (extracted_xlims[2] + 0.01), y = y, 
                                    label=Country, 
                                    color = Region), 
                               size=2,
                               hjust = "left")


ggt2_hm <- ggt2_hm + geom_point(data = ggt2_hm$data[which(!is.na(ggt2_hm$data$novel_seq)),],
                          aes( x = (extracted_xlims[2] - 0.13), y = y), 
                          color = "black",
                          size=1, shape = 8)



ggt2_hm <- ggt2_hm + xlim((extracted_xlims[1]), (extracted_xlims[2]+0.1))



#### Coloring of mosquito genus and region
ggt2_hm <- ggt2_hm + scale_fill_manual(breaks=c('Anopheles',
                                                'Armigeres',
                                                'Coquillettidia',
                                                'Culiseta',
                                                'Culex',
                                                'Aedes',
                                                
                                                'Culicidae',
                                                'Arthropoda',
                                                
                                                'na',
     
                                                "North America and Caribbean", 
                                                "South America", 
                                                "Europe", 
                                                "West Africa",
                                                "Southeast Africa",
                                                "Southwest Asia",
                                                "Southeast Asia"                                      
                                                ),
                                       values=c("black",
                                                "#88CCEE",
                                                "#117733",
                                                "#DDCC77",
                                                
                                                "#332288",
                                                
                                                "#882255",
                                                
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

ggt2_hm <- ggt2_hm + scale_color_manual(breaks=c('Anopheles',
                                                 'Armigeres',
                                                 'Coquillettidia',
                                                 'Culiseta',
                                                 'Culex',
                                                 'Aedes',
                                                 
                                                 'Culicidae',
                                                 'Arthropoda',
                                                 
                                                 'na',
                                                 
                                                 "North America and Caribbean", 
                                                 "South America", 
                                                 "Europe", 
                                                 "West Africa",
                                                 "Southeast Africa",
                                                 "Southwest Asia",
                                                 "Southeast Asia"                                      
),
values=c("black",
         "#88CCEE",
         "#117733",
         "#DDCC77",
         
         "#332288",
         
         "#882255",
         
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
ggt2_hm <- ggt2_hm + ylim((extracted_ylims[1]-5), (extracted_ylims[2]+5))
extracted_xlims1 <- layer_scales(ggt2_hm)$x$range$range
ggt2_hm <- ggt2_hm + xlim((extracted_xlims1[1]), (extracted_xlims1[2]+0.1))
ggt2_hm <- ggt2_hm + theme_tree(legend.position="none")
ggt2_hm <- ggt2_hm + theme(plot.margin = unit(c(0, 0, 0, 0), "inch"))

ggt2_hm_legend <- ggt2_hm + theme_tree(legend.position="bottom")
# 
# print(ggt2_hm)
ggt2_hm_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/Sobeli_clade_aa_v2311_dedup_aln.ggtree_hm.pdf"
ggsave(plot = ggt2_hm, filename = ggt2_hm_file, device = "pdf", width = 10, height = 14, units = "cm")

# ggt2_hm_file_legend <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/Sobeli_clade_aa_v2311_dedup_aln.ggtree_hm_legend.pdf"
# ggsave(plot = ggt2_hm_legend, filename = ggt2_hm_file_legend, device = "pdf", width = 50, height = 50, units = "cm")














####################### 
####################### TREES. SUPPLEMENTARY FIGURE
####################### 
tip_mdf_brief$nt_seq_ID
tip_mdf_brief$var <- tip_mdf_brief$host_taxon1

genus_matrix <- tip_mdf_brief %>% select(nt_seq_ID, host_taxon1, var) %>% 
  spread(key = host_taxon1, value = var, fill = NA) %>% 
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
ggt2_hm <- ggt2_hm + scale_fill_manual(breaks=c('Anopheles',
                                                'Armigeres',
                                                'Coquillettidia',
                                                'Culiseta',
                                                'Culex',
                                                'Aedes',
                                                
                                                'Culicidae',
                                                'Arthropoda',
                                                
                                                'na',
                                                
                                                "North America and Caribbean", 
                                                "South America", 
                                                "Europe", 
                                                "West Africa",
                                                "Southeast Africa",
                                                "Southwest Asia",
                                                "Southeast Asia"                                      
),
values=c("black",
         "#88CCEE",
         "#117733",
         "#DDCC77",
         
         "#332288",
         
         "#882255",
         
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
ggt_sup2_hm_file <- "/full_path_to/wd/RdRp_scan/analysis/phylo/Sobeli_clade/Sobeli_clade_aa_v2311_dedup_aln.ggtree_hm_supplementary.pdf"
ggsave(plot = ggt_sup2_hm, filename = ggt_sup2_hm_file, device = "pdf", width = 21, height = 60, units = "cm")
