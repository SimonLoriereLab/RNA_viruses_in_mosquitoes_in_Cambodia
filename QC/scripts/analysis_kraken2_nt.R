library(tidyverse)


########### Run this part on MAESTRO HPC


# path_to_inout <- "/full_path_to/wd/QC/analysis/Kraken2/"
# date <- str_remove_all(Sys.Date(), "-")
# 
# 
# k2df1 <- read_delim(paste0(path_to_inout, "kraken2_reports_CMVM_20211116.txt.gz"), delim = "\t",
#                     num_threads = 95) ##### this file is 13 Gb - not in github repo
# 
# k2df1_short <- k2df1[which(k2df1$taxid %in% c(0, ### Unclassified
#                                               1, ### Root
#                                               2157, ### Archaea
#                                               2,	### Bacteria
#                                               2759, ### Eukaryotes
#                                               7157, ### Culicidae
#                                               10239 ### Viruses
#                                               )),]
# 
# write.table(x = k2df1_short, file = gzfile(paste0(path_to_inout, "kraken2_reports_CMVM_20220907.txt.gz_short.txt.gz")),
#             append = F, quote = F, sep = "\t", row.names = F, col.names = T)






########### Run this part on desktop

path_to_inout <- "/full_path_to/wd/QC/analysis/Kraken2/"
date <- str_remove_all(Sys.Date(), "-")
k2df1_short <- read_delim(paste0(path_to_inout, "kraken2_reports_CMVM_20220907.txt.gz_short.txt.gz"), delim = "\t",
                          num_threads = 95)
mdf <- read_delim("/full_path_to/wd/QC/metadata/CMVM_full_metadata_updated_20211021.tsv", delim = "\t") %>%
  select(sample_id = Sample, Batch, LibPool, NTC, Location, genus_species)


tdf <- read_delim("/full_path_to/wd/QC/analysis/Trimmomatic/trimmomatic_summary_CMVM_20211122.txt", delim = "\t")

tdf1 <- tdf[which(tdf$flag_poor_read_qual == 0),]


length(unique(tdf$sample_id))


#### prepare datasets
k2df1_short <- left_join(k2df1_short, mdf, by = "sample_id")
k2df1_short$percent_frags <- as.numeric(str_trim(k2df1_short$percent_frags, side = "both"))
k2df1_short1 <- k2df1_short

#### Removing Bacteria, Archaea, Eukaryota, Viruses from the root
k2df1_short1$taxon_name <- str_trim(k2df1_short1$taxon_name, side = "both")
k2df1_short1[which(k2df1_short1$taxid == "1"), ]$percent_frags <- (k2df1_short[which(k2df1_short$taxid == "1"), ]$percent_frags - 
                                                                     k2df1_short[which(k2df1_short$taxid == "2157"), ]$percent_frags -
                                                                     k2df1_short[which(k2df1_short$taxid == "2"), ]$percent_frags - 
                                                                     k2df1_short[which(k2df1_short$taxid == "2759"), ]$percent_frags - 
                                                                     k2df1_short[which(k2df1_short$taxid == "10239"), ]$percent_frags)

k2df1_short1[which(k2df1_short1$taxid == "1"), ]$num_frags <- (k2df1_short[which(k2df1_short$taxid == "1"), ]$num_frags - 
                                                                     k2df1_short[which(k2df1_short$taxid == "2157"), ]$num_frags -
                                                                     k2df1_short[which(k2df1_short$taxid == "2"), ]$num_frags - 
                                                                     k2df1_short[which(k2df1_short$taxid == "2759"), ]$num_frags - 
                                                                     k2df1_short[which(k2df1_short$taxid == "10239"), ]$num_frags)
k2df1_short1[which(k2df1_short1$taxid == "1"), ]$taxon_name <- "Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)"

### excluding Culicidae from Eukaryotes
k2df1_short1[which(k2df1_short1$taxid == "2759"), ]$percent_frags <- (k2df1_short[which(k2df1_short$taxid == "2759"), ]$percent_frags - k2df1_short[which(k2df1_short$taxid == "7157"), ]$percent_frags)
k2df1_short1[which(k2df1_short1$taxid == "2759"), ]$num_frags <- (k2df1_short[which(k2df1_short$taxid == "2759"), ]$num_frags - k2df1_short[which(k2df1_short$taxid == "7157"), ]$num_frags)
k2df1_short1$taxon_name <- str_trim(k2df1_short1$taxon_name, side = "both")
k2df1_short1[which(k2df1_short1$taxid == "2759"), ]$taxon_name <- "Eukaryota (-Culicidae)"








##
k2df1_short2 <- k2df1_short1 %>% select(num_frags, taxon_name) %>% group_by(taxon_name) %>% summarize(tax_num_frags = sum(num_frags))
sum(k2df1_short1$num_frags) == sum(k2df1_short2$tax_num_frags) ### check
k2df1_short2$perc <- 100*k2df1_short2$tax_num_frags/sum(k2df1_short2$tax_num_frags)
#####
k2df1_short2$taxon_name <- factor(k2df1_short2$taxon_name,
                                  levels = c('Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)','Archaea','Viruses',
                                             'Bacteria','unclassified', 'Culicidae', 'Eukaryota (-Culicidae)'),
                                  ordered = T)
cat(paste0(k2df1_short2$taxon_name, collapse = "', '"))


k2df1_short2$label <- paste0(round(k2df1_short2$tax_num_frags/(10^6)),
                             "M - ",
                             round(k2df1_short2$perc, digits = 1), "%")
k2df1_short2[which(k2df1_short2$taxon_name %in% c('Culicidae', 'unclassified')),]$label <- paste0(round(k2df1_short2[which(k2df1_short2$taxon_name %in% c('Culicidae', 'unclassified')),]$tax_num_frags/(10^9), 
                                                                                                        digits = 1),
                             "B - ",
                             round(k2df1_short2[which(k2df1_short2$taxon_name %in% c('Culicidae', 'unclassified')),]$perc, digits = 1), "%")

k2df1_short2[which(k2df1_short2$taxon_name %in% c('Archaea')),]$label <- paste0(round(k2df1_short2[which(k2df1_short2$taxon_name %in% c('Archaea')),]$tax_num_frags/(10^3)),
                                                                                                  "K - ",
                                                                                                  round(k2df1_short2[which(k2df1_short2$taxon_name %in% c('Archaea')),]$perc, digits = 1), "%")



p1 <- ggplot(k2df1_short2, aes(x = "", y = perc, fill = taxon_name)) +
  geom_col(color = "black") +
  geom_label(aes(label = label, color  = taxon_name),
             size = 2,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +

  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c('Viruses'='#a50f15', 'Culicidae' = '#de2d26', 'Eukaryota (-Culicidae)' = '#fb6a4a',
                               'Archaea' = '#fcae91', 'Bacteria' = '#fee5d9','Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)' = 'grey30',
                               'unclassified' = 'grey80')) +
  scale_color_manual(values = c('Viruses'='white', 'Culicidae' = 'white', 'Eukaryota (-Culicidae)' = 'white',
                                'Archaea' = 'black', 'Bacteria' = 'black', 'Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)' = 'white',
                               'unclassified' = 'black')) +
  coord_polar(theta = "y") +
  theme_void(base_size = 7)

p1

ggsave(filename = paste0(path_to_inout, "pie_", date, ".pdf"),
       plot = p1, device = "pdf", width = 10, height = 4, units = "cm")









################################################################################
##### Supplementary pies for Mosquitoes and NTCs

##
k2df1_short3 <- k2df1_short1 %>% select(num_frags, taxon_name, NTC) %>% group_by(taxon_name, NTC) %>% summarize(tax_num_frags = sum(num_frags))
k2df1_short3$perc <- NA
k2df1_short3[which(k2df1_short3$NTC == 0), ]$perc <- 100*k2df1_short3[which(k2df1_short3$NTC == 0), ]$tax_num_frags/sum(k2df1_short3[which(k2df1_short3$NTC == 0), ]$tax_num_frags)
k2df1_short3[which(k2df1_short3$NTC == 1), ]$perc <- 100*k2df1_short3[which(k2df1_short3$NTC == 1), ]$tax_num_frags/sum(k2df1_short3[which(k2df1_short3$NTC == 1), ]$tax_num_frags)
#####
k2df1_short3$taxon_name <- factor(k2df1_short3$taxon_name,
                                  levels = c('Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)','Archaea','Viruses','Bacteria','unclassified', 'Culicidae', 'Eukaryota (-Culicidae)'),
                                  ordered = T)



k2df1_short3$label <- paste0(round(k2df1_short3$tax_num_frags/(10^6)),
                             "M - ",
                             round(k2df1_short3$perc, digits = 1), "%")
k2df1_short3[which(k2df1_short3$taxon_name %in% c('Culicidae', 'unclassified')),]$label <- paste0(round(k2df1_short3[which(k2df1_short3$taxon_name %in% c('Culicidae', 'unclassified')),]$tax_num_frags/(10^9),
                                                                                                        digits = 1),
                                                                                                  "B - ",
                                                                                                  round(k2df1_short3[which(k2df1_short3$taxon_name %in% c('Culicidae', 'unclassified')),]$perc, digits = 1), "%")

k2df1_short3[which(k2df1_short3$taxon_name %in% c('Archaea')),]$label <- paste0(round(k2df1_short3[which(k2df1_short3$taxon_name %in% c('Archaea')),]$tax_num_frags/(10^3)),
                                                                                                  "K - ",
                                                                                                  round(k2df1_short3[which(k2df1_short3$taxon_name %in% c('Archaea')),]$perc, digits = 1), "%")


k2df1_short3[which(k2df1_short3$taxon_name %in%
                     c('Archaea','Viruses', 'Culicidae') &
                     k2df1_short3$NTC == 1),]$label <- paste0(round(k2df1_short3[which(k2df1_short3$taxon_name %in%
                                                                                         c('Archaea','Viruses', 'Culicidae') &
                                                                                         k2df1_short3$NTC == 1),]$tax_num_frags/(10^3)),
                                                                                                  "K - ",
                                                                                                  round(k2df1_short3[which(k2df1_short3$taxon_name %in%
                                                                                                                             c('Archaea','Viruses', 'Culicidae') &
                                                                                                                             k2df1_short3$NTC == 1),]$perc, digits = 1), "%")
k2df1_short3[which(k2df1_short3$taxon_name %in%
                     c('unclassified') &
                     k2df1_short3$NTC == 1),]$label <- paste0(round(k2df1_short3[which(k2df1_short3$taxon_name %in%
                                                                                         c('unclassified') &
                                                                                         k2df1_short3$NTC == 1),]$tax_num_frags/(10^6)),
                                                              "M - ",
                                                              round(k2df1_short3[which(k2df1_short3$taxon_name %in%
                                                                                         c('unclassified') &
                                                                                         k2df1_short3$NTC == 1),]$perc, digits = 1), "%")

k2df1_short3_m <- k2df1_short3[which(k2df1_short3$NTC == 0),]

p2 <- ggplot(k2df1_short3_m, aes(x = "", y = perc, fill = taxon_name)) +
  geom_col(color = "black") +
  geom_label(aes(label = label, color  = taxon_name),
             size = 2,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +

  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c('Viruses'='#a50f15', 'Culicidae' = '#de2d26', 'Eukaryota (-Culicidae)' = '#fb6a4a',
                               'Archaea' = '#fcae91', 'Bacteria' = '#fee5d9','Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)' = 'grey30',
                               'unclassified' = 'grey80')) +
  scale_color_manual(values = c('Viruses'='white', 'Culicidae' = 'white', 'Eukaryota (-Culicidae)' = 'white',
                                'Archaea' = 'black', 'Bacteria' = 'black', 'Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)' = 'white',
                                'unclassified' = 'black')) +
  coord_polar(theta = "y") +
  theme_void(base_size = 7)

p2

ggsave(filename = paste0(path_to_inout, "pie_mosq_", date, ".pdf"),
       plot = p2, device = "pdf", width = 15, height = 7, units = "cm")



k2df1_short3_n <- k2df1_short3[which(k2df1_short3$NTC == 1),]


p3 <- ggplot(k2df1_short3_n, aes(x = "", y = perc, fill = taxon_name)) +
  geom_col(color = "black") +
  geom_label(aes(label = label, color  = taxon_name),
             size = 2,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +

  guides(fill = guide_legend(title = "")) +
  scale_fill_manual(values = c('Viruses'='#a50f15', 'Culicidae' = '#de2d26', 'Eukaryota (-Culicidae)' = '#fb6a4a',
                               'Archaea' = '#fcae91', 'Bacteria' = '#fee5d9','Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)' = 'grey30',
                               'unclassified' = 'grey80')) +
  scale_color_manual(values = c('Viruses'='white', 'Culicidae' = 'white', 'Eukaryota (-Culicidae)' = 'white',
                                'Archaea' = 'black', 'Bacteria' = 'black', 'Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)' = 'white',
                                'unclassified' = 'black')) +
  coord_polar(theta = "y") +
  theme_void(base_size = 7)

p3

ggsave(filename = paste0(path_to_inout, "pie_NTC_", date, ".pdf"),
       plot = p3, device = "pdf", width = 15, height = 7, units = "cm")



##### Supplementary stacked proportions for individual samples
k2df1_short1$Batch <- floor(k2df1_short1$Batch)

k2df1_short1$sample_id <- factor(k2df1_short1$sample_id,
                                  levels = sort(unique(k2df1_short1$sample_id)),
                                  ordered = T)

k2df1_short1$taxon_name <- factor(k2df1_short1$taxon_name,
                                  levels = c('Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)','unclassified',
                                             'Bacteria', 'Archaea',
                                             'Eukaryota (-Culicidae)',
                                             'Culicidae', 'Viruses'),
                                  ordered = T)
k2df1_short1$alpha <- "m"
k2df1_short1[which(k2df1_short1$NTC == 1),]$alpha <- "n"


unique(k2df1_short1[which(k2df1_short1$Batch == 17),]$sample_id)


p4 <-  ggplot()
p4 <- p4 + geom_col(data = k2df1_short1, mapping = aes(x = sample_id, y = percent_frags,
                                                       fill = taxon_name, color = alpha), size = 0.25)
p4 <- p4 + scale_fill_manual(values = c('Viruses'='#a50f15', 'Culicidae' = '#de2d26', 'Eukaryota (-Culicidae)' = '#fb6a4a',
                             'Archaea' = '#fcae91', 'Bacteria' = '#fee5d9','Root (-Bacteria, -Archaea, -Eukaryota, -Viruses)' = 'grey30',
                             'unclassified' = 'grey80'))
p4 <- p4 + scale_color_manual(values = c("white","black"), guide  = "none")
# p4 <- p4 + scale_alpha_manual(values = c(1,0.5), guide  = "none")
p4 <- p4 + theme_classic(base_size = 5)
p4 <- p4 + facet_wrap(facets = "Batch", ncol = 4, scales = "free_x")
p4 <- p4 + theme(legend.position = "none",
                 axis.text.x = element_text(angle = 90))
# p4

ggsave(filename = paste0(path_to_inout, "stacked_bar_all_", date, ".pdf"),
       plot = p4, device = "pdf", width = 20, height = 20, units = "cm")

