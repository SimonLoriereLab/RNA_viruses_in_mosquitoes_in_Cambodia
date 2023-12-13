getwd()
setwd("/full_path_to/wd/before_NGS/")

library(tidyverse)
library(readxl)
library(viridis)
library(measurements)
library(ggrepel)
library(gridExtra)
library(cowplot)

mdf3 <- read.delim("./metadata/FSPI_samples_Cambodia_mosquito_viral_menagenomics.txt")


##### plot on the map
dir.create("plots")
cambodia <- raster::getData("GADM",country="Cambodia",level=1)
camb <- fortify(cambodia)
cambo <- cbind(camb, cambodia@data[camb$id,])
province <- c("other divisions","Batdâmbâng","Kâmpóng Thum","Krong Pailin","Krong Preah Sihanouk","Pouthisat","Preah Vihéar")
province_num <- c("12", "4", "5", "9", "22", "10")
province_palette <- c('grey90', '#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc')
cambo$division <- ""
cambo[which(!cambo$id %in% province_num),]$division <- "other divisions"
cambo[which(cambo$id %in% province_num),]$division <- cambo[which(cambo$id %in% province_num),]$NAME_1 
unique(cambo$division)

cambo$division <- as.factor(cambo$division)

map3 <- ggplot(data = mdf3)
map3 <- map3 + geom_polygon(data = cambo, aes(x = long, y = lat, 
                                              group = group, 
                                              fill = division),colour = "white")
map3 <- map3 + coord_fixed(ratio=1)
map3 <- map3 + geom_point(data = mdf3, aes(x = as.numeric(long), y = as.numeric(lat)), size = 1)
map3 <- map3 + theme_minimal(base_size = 9, base_family = "Helvetica")
map3 <- map3 + scale_fill_manual(breaks = province, 
                                 values = c('grey90', '#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc'))
map3 <- map3 + theme(legend.position = "bottom")
ggsave(filename = "plots/map_Cambodia_samples.pdf", plot = map3, 
       device = "pdf",width = 20, height = 30, units = "cm", dpi = 300)


### plot on the map by species
species_list <- unique(mdf3$genus_species)
unique(mdf3$trip)

summary_mdf3 <- mdf3 %>% 
  group_by(lat, long, genus_species, site) %>%
  summarise(n = n())

dir.create("plots/samles_on_the_map_by_species_counts")
for (s in species_list){
  map3 <- ggplot(data = summary_mdf3[which(summary_mdf3$genus_species == s),])
  map3 <- map3 + geom_polygon(data = cambo, aes(x = long, y = lat, 
                                                group = group, 
                                                fill = division),colour = "white")
  map3 <- map3 + coord_fixed(ratio=1)
  map3 <- map3 + geom_point(data = summary_mdf3[which(summary_mdf3$genus_species == s),],
                            aes(x = as.numeric(long),y = as.numeric(lat)), 
                            size = 2, shape = 10, color = "firebrick2")
  map3 <- map3 + geom_text_repel(data = summary_mdf3[which(summary_mdf3$genus_species == s),],
                                 mapping = aes(x = as.numeric(long), y = as.numeric(lat), label = n),
                                 seed = 7)
  map3 <- map3 + theme_minimal(base_size = 9, base_family = "Helvetica")
  map3 <- map3 + ggtitle(label = s)
  map3 <- map3 + scale_fill_manual(breaks = province, 
                                   values = c('grey90', '#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc'))
  map3 <- map3 + theme(legend.position = "bottom")
  ggsave(filename = paste0("plots/samles_on_the_map_by_species_counts/", s, ".pdf"), plot = map3, 
         device = "pdf", width = 20, height = 28, units = "cm", dpi = 300)
}


summary_mdf3_sites <- mdf3 %>% 
  group_by(lat, long, site) %>%
  summarise(n = n())

map4 <- ggplot(data = summary_mdf3_sites)
map4 <- map4 + geom_polygon(data = cambo, aes(x = long, y = lat, 
                                              group = group, 
                                              fill = division),colour = "white")
map4 <- map4 + coord_fixed(ratio=1)
map4 <- map4 + geom_point(data = summary_mdf3_sites, aes(x = as.numeric(long), y = as.numeric(lat)), size = 1)
map4 <- map4 + geom_text_repel(data = summary_mdf3_sites, 
                               mapping = aes(x = as.numeric(long), y = as.numeric(lat), label = site), seed = 7,  max.overlaps = Inf)
map4 <- map4 + theme_minimal(base_size = 9, base_family = "Helvetica")
map4 <- map4 + scale_fill_manual(breaks = province, 
                                 values = c('grey90', '#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc'))
map4 <- map4 + theme(legend.position = "bottom")
ggsave(filename = "plots/map_Cambodia_sites.pdf", plot = map4, 
       device = "pdf",width = 20, height = 30, units = "cm", dpi = 300)


###### Select species that would be a priority for sequencing based on the overall representation across sites and regions
selected_sp <- c("Aedes aegypti","Aedes albopictus","Culex quinquefasciatus",
  "Mansonia indiana","Culex sinensis","Culex vishnui.g",
  "Culex bitaeniorhynchus","Coquillettidia crassipes",
  "Armigeres subalbatus","Culex gelidus")

mdf3$seq_priority <- "low"
mdf3[which(mdf3$genus_species %in% selected_sp),]$seq_priority <- "high"
setdiff(selected_sp, unique(mdf3[which(mdf3$seq_priority == "high"),]$genus_species))
write.table(x = mdf3, file = "metadata/full_df_seq_priority_annotated.txt", 
                        append = FALSE, quote = FALSE, sep = "\t'", row.names = FALSE, col.names = TRUE)

### import into a Google sheet for inspection, download and reimport below
####### update the first batch samples used
mdf6 <- read.delim("metadata/CMVM_Mosquito_viromics_cambodia - priority annotated full sample list.tsv", header = T)
bldf <- read.delim("metadata/CMVM_Mosquito_viromics_cambodia - batch_lib.tsv", header = T)
curdf <- full_join(bldf, mdf6, by = "Sample", .keep_all = T)

write.table(x = curdf, file = "full_df_seq_priority_annotated_batch1_lp1.txt",
            append = FALSE, quote = FALSE, sep = "\t'",row.names = FALSE, col.names = TRUE)

####### prepare list of all further batches and libpools

length(curdf[which(curdf$seq_priority == "high"),]$Sample)

hpdf_start_b2 <- curdf[which(curdf$seq_priority == "high" & is.na(curdf$Batch)),]

##### DO NOT rerun the randomization of the samples code below
##### hpdf_start_b2_rand1 <- hpdf_start_b2[sample(nrow(hpdf_start_b2)),]    ##### did this only once

# length(hpdf_start_b2_rand1$Batch)%/%30
# length(hpdf_start_b2_rand1$Batch)%%30
# 
# length(c(rep(2:16, each=30), rep(17, 8)))
# hpdf_start_b2_rand1$Batch <- paste0(as.character(c(rep(2:16, each=30), rep(17, 8))),".", as.character(c(rep(rep(1:2, each=15), 15), rep(1, 8))))

#order by tube by box and obviously by batch for convenience

# hpdf_start_b2_rand11 <- hpdf_start_b2_rand1[
#   with(hpdf_start_b2_rand1, order(as.numeric(Batch), Location, as.numeric(new_box), as.numeric(new_eppendorf))),
#   ]
# 
# 
# hpdf_start_lp2_first7 <- curdf[which(curdf$seq_priority == "high" & is.na(curdf$Batch)==FALSE & is.na(curdf$LibPool)),]
# 
# hpdf_start_lp2_first7$Batch <- as.character(hpdf_start_lp2_first7$Batch)
# 
# hpdf_start_lp2 <- bind_rows(hpdf_start_lp2_first7, hpdf_start_b2_rand11)
# 
# hpdf_start_lp2 
# 
# length(hpdf_start_lp2$Batch)%/%22
# length(hpdf_start_lp2$Batch)%%22
# 
# hpdf_start_lp2$LibPool <- as.character(c(rep(2:22, each=22), rep(23, 3)))
# 
# 
# hpdf_lp1 <- curdf[which(curdf$seq_priority == "high" & curdf$LibPool == 1),]
# 
# hpdf_lp1$Batch <- as.character(hpdf_lp1$Batch)
# hpdf_lp1$LibPool <- as.character(hpdf_lp1$LibPool)
# 
# hpdf_all <- bind_rows(hpdf_lp1, hpdf_start_lp2)
# 
# lowdf <- curdf[which(curdf$seq_priority == "low" ),]
# 
# lowdf$Batch <- "0.0"
# lowdf$LibPool <- "0"
# 
# hpdf_all_and_low <- bind_rows(hpdf_all, lowdf)
# 
# setdiff(curdf$Sample, hpdf_all_and_low$Sample)
# 
# write.table(x = hpdf_all_and_low, file = "metadata/all_samples_batches_libPools_assigned.txt",
#             append = F, quote = F, sep = "\t", row.names = F, col.names = T)



##### import df generated above
hpdf_all_and_low <- read_delim("metadata/all_samples_batches_libPools_assigned.txt", delim = "\t")
hpdf_all <- hpdf_all_and_low[which(hpdf_all_and_low$seq_priority == "high"),]

checkdf1 <- hpdf_all %>% select(Batch, genus_species) %>%
  group_by(Batch, genus_species) %>%
  summarise(total = n())

p1 <- ggplot(checkdf1, aes(x = genus_species, y = Batch, fill = total))+
  geom_tile()+
  geom_text(aes(label = total))+
  scale_fill_distiller(palette = "RdPu", direction = 1)+
  theme_minimal(base_size = 9)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

checkdf2 <- hpdf_all %>% select(LibPool, genus_species) %>%
  group_by(LibPool, genus_species) %>%
  summarise(total = n())

p2 <- ggplot(checkdf2, aes(x = genus_species, y = LibPool, fill = total))+
  geom_tile()+
  geom_text(aes(label = total))+
  scale_fill_distiller(palette = "RdPu", direction = 1)+
  theme_minimal(base_size = 9)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

checkdf3 <- hpdf_all %>% select(Batch, Location) %>%
  group_by(Batch, Location) %>%
  summarise(total = n())

p3 <- ggplot(checkdf3, aes(x = Location, y = Batch, fill = total))+
  geom_tile()+
  geom_text(aes(label = total))+
  scale_fill_distiller(palette = "RdPu", direction = 1)+
  theme_minimal(base_size = 9)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

checkdf4 <- hpdf_all %>% select(LibPool, Location) %>%
  group_by(LibPool, Location) %>%
  summarise(total = n())

p4 <- ggplot(checkdf4, aes(x = Location, y = LibPool, fill = total))+
  geom_tile()+
  geom_text(aes(label = total))+
  scale_fill_distiller(palette = "RdPu", direction = 1)+
  theme_minimal(base_size = 9)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p5 <- varcovplot_combined <- plot_grid(p1, p3, p2, p4, nrow = 2, ncol = 2)
save_plot(paste0("plots/check_sample_randomization.eps"), p5,  nrow = 2, ncol = 2,  base_height = 7, base_width = 7)


