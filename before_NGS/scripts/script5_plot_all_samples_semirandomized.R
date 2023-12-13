getwd()
setwd("/full_path_to/wd/before_NGS/")
library(tidyverse)


mdf <- read_delim("metadata/CMVM_full_metadata_updated_20211021.tsv", delim = '\t')

mdf1 <- mdf[which(!is.na(mdf$LibPool)), ]

mdf1$Batch <- factor(as.character(mdf1$Batch), ordered = T, levels = as.character(sort(unique(mdf1$Batch))))

checkdf3 <- mdf1 %>% select(Batch, Location) %>%
  group_by(Batch, Location) %>%
  summarise(total = n())

p3 <- ggplot(checkdf3, aes(x = Location, y = Batch, fill = total))+
  geom_tile()+
  geom_text(aes(label = total))+
  scale_fill_distiller(palette = "RdPu", direction = 1)+
  theme_minimal(base_size = 9)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p3)
ggsave(plot = p3, filename = "plots/Batch_x_Loc.pdf", device = "pdf", width = 15, height = 15, units = "cm")


checkdf4 <- mdf1 %>% select(genus_species, Batch) %>%
  group_by(Batch, genus_species) %>%
  summarise(total = n())

p4 <- ggplot(checkdf4, aes(x = genus_species, y = Batch, fill = total))+
  geom_tile()+
  geom_text(aes(label = total))+
  scale_fill_distiller(palette = "RdPu", direction = 1)+
  theme_minimal(base_size = 9)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p4)
ggsave(plot = p4, filename = "plots/Batch_x_sp.pdf", device = "pdf", width = 15, height = 15, units = "cm")

