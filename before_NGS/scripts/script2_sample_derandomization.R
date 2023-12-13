getwd()
setwd("/full_path_to/wd/before_NGS/")

library(tidyverse)
library(readxl)
library(viridis)
library(measurements)
library(ggrepel)
library(gridExtra)
library(cowplot)


all_df <- read.delim2("metadata/CMVM_Mosquito_viromics_cambodia_batch5 - priority annotated full sample list_05102020.tsv", sep = "\t")

lpdf <- all_df[which(all_df$seq_priority == "low"),]
processed_hpdf <- all_df[which(all_df$LibPool < 6 & all_df$seq_priority == "high"),]
lp17_hpdf <- all_df[which(all_df$LibPool == 17 & all_df$seq_priority == "high"),] ### special case, reason related to indexing kits

# here samples that have not been processed. 
# I am reordering and changing batches to process them sequentially
# according to the tube ordering to retrieve samples faster to prevent RNA degradation

todo_hpdf <- all_df[which(all_df$LibPool > 5 & all_df$LibPool < 17 & all_df$seq_priority == "high"),]


todo_hpdf1 <- todo_hpdf[order(todo_hpdf$Location, todo_hpdf$new_eppendorf),]

length(todo_hpdf1$Batch)%/%30
length(todo_hpdf1$Batch)%%30

length(c(rep(6:16, each=30)))

todo_hpdf1$Batch <- paste0(as.character(c(rep(6:16, each=30))),".", as.character(c(rep(rep(1:2, each=15), 11))))

todo_hpdf1$LibPool <- str_remove(todo_hpdf1$Batch, "\\.[[:digit:]]")

todo_hpdf1$LibPool <- as.numeric(todo_hpdf1$LibPool)
todo_hpdf1$Batch <- as.character(todo_hpdf1$Batch)

hpdf1 <- bind_rows(processed_hpdf, todo_hpdf1, lp17_hpdf, lpdf)

hpdf_all <- hpdf1
hpdf_all$Batch <- as.numeric(hpdf_all$Batch)

# write.table(x = hpdf_all, file = "metadata/all_samples_batches_libPools_reassigned_05102020_tick.txt",
#             append = F, quote = F, sep = "\t'", row.names = F, col.names = T)

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

p5 <- varcovplot_combined <- plot_grid(p2, p4, nrow = 1, ncol = 2)
save_plot(paste0("plots/check_sample_randomization_back_to_sequential.eps"), p5,  nrow = 1, ncol = 2,  base_height = 7, base_width = 7)




