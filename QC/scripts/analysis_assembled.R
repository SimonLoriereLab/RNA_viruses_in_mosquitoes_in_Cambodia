library(tidyverse)
library(ggrepel)
library(ggridges)
library(GGally)
library(viridisLite)


path_to_inout <- "/full_path_to/wd/QC/analysis/metaspades/"
date <- str_remove_all(Sys.Date(), "-")

# #### preprocessing
### DONE
# mdf <- read_delim("/full_path_to/wd/QC/metadata/CMVM - full_metadata_20211006.tsv", delim = "\t") %>%
#   select(sample_id = Sample, Batch, LibPool, NTC, Location, genus_species)
# 
# df_denovo <- read_delim(file = paste0(path_to_inout, "metaspades_scraped_stats_CMVM_20211013.txt.gz"), delim = "\t")
# 
# cat(paste0(colnames(df_denovo), collapse = ", "))
# 
# df_denovo <- left_join(df_denovo, mdf, by = "sample_id")
# 
# df_trim <- read_delim(file = "/full_path_to/wd/QC/analysis/Trimmomatic/trimmomatic_scraped_stats_CMVM_20211005.txt",
#                       delim = "\t") %>% select(sample_id, total_num_of_PE_reads, num_qualtrim_kept_PE_reads)
# 
# df_denovo <- left_join(df_denovo, df_trim, by = "sample_id") %>% distinct()

# write.table(x = df_denovo, file = gzfile(paste0(path_to_inout, "metaspades_scraped_stats_CMVM_20211013_with_trim_and_meta.txt.gz")),
#             append = F, quote = F, sep = "\t", row.names = F, col.names = T)


df_denovo <- read_delim(paste0(path_to_inout, "metaspades_scraped_stats_CMVM_20211013_with_trim_and_meta.txt.gz"), delim = "\t") %>% distinct()

df_denovo$Batch <-  factor(as.character(df_denovo$Batch), levels = as.character(sort(unique(df_denovo$Batch))))
class(df_denovo$Batch)
levels(df_denovo$Batch)

df_denovo$LibPool <-  factor(as.character(df_denovo$LibPool), levels = as.character(sort(unique(df_denovo$LibPool))))
class(df_denovo$LibPool)
levels(df_denovo$LibPool)




df_denovo$NTC <- as.factor(as.character(df_denovo$NTC))



# These are just very low coverage in 33 contigs from 31 samples from 19 batches. Mask them from plots on log scale.
sum(df_denovo$contig_kmer_cov == 0)
df_denovo_zero_cov <- df_denovo[which(df_denovo$contig_kmer_cov == 0),]
sort(unique(df_denovo_zero_cov$sample_id))
length(sort(unique(df_denovo_zero_cov$Batch)))

df_denovo_non_zero <-  df_denovo[which(df_denovo$contig_kmer_cov > 0),]

#### plots


# kmer covs 
### DONE
p0 <- ggplot()
p0 <- p0 + geom_histogram(data = df_denovo_non_zero, mapping = aes(x = log10(contig_kmer_cov)), bins = 100)
p0 <- p0 + theme_minimal(base_size = 11)
#p0

ggsave(filename = paste0(path_to_inout, "kmer_cov_", date, ".pdf"),
       plot = p0, device = "pdf", width = 15, height = 10, units = "cm")


### ridge plot

df_denovo_nz_summary <- df_denovo_non_zero %>% 
  group_by(sample_id, total_num_of_PE_reads, num_qualtrim_kept_PE_reads,
           Batch, LibPool, NTC, Location, genus_species) %>% 
  summarise(average_contig_len = mean(contig_len),
            average_contig_kmer_cov = mean(contig_kmer_cov),
            total_num_of_contigs = n(), lensum_of_contigs= sum(contig_len)) %>%
  arrange(desc(average_contig_len))

df_denovo_non_zero <- left_join(df_denovo_non_zero, df_denovo_nz_summary)


df_denovo_non_zero$sample_id_len_ord <- factor(as.character(df_denovo_non_zero$sample_id), 
                                               ordered = T, 
                                               levels = df_denovo_nz_summary$sample_id)


df_denovo_nz_summary1 <- df_denovo_nz_summary %>%
  arrange(desc(average_contig_kmer_cov))



df_denovo_non_zero$sample_id_kmer_cov_ord <- factor(as.character(df_denovo_non_zero$sample_id), 
                                               ordered = T, 
                                               levels = df_denovo_nz_summary1$sample_id)

df_denovo_nz_summary$prop_qualtrim_kept_PE_reads <- df_denovo_nz_summary$num_qualtrim_kept_PE_reads/df_denovo_nz_summary$total_num_of_PE_reads


p3 <- ggplot(data = df_denovo_non_zero, aes(x = log10(contig_len), y = sample_id_len_ord))
p3 <- p3 + geom_density_ridges(aes(fill = LibPool), scale = 2, size = 0.5, stat = "binline", bins = 100)
p3 <- p3 + theme_ridges(font_size = 6, grid = FALSE)
#p3 <- p3 + facet_grid(genus_species~., scales = "free_y", space = "free_y")
ggsave(filename = paste0(path_to_inout, "ridge_cov_len_", date, ".pdf"),
        plot = p3, device = "pdf", width = 10, height = 100, units = "cm")

p4 <- ggplot(data = df_denovo_non_zero, aes(x = log10(contig_kmer_cov), y = sample_id_len_ord))
p4 <- p4 + geom_density_ridges(aes(fill = LibPool), scale = 2, size = 0.5, stat = "binline", bins = 100)
p4 <- p4 + theme_ridges(font_size = 6, grid = FALSE)
#p4 <- p4 + facet_grid(genus_species~., scales = "free_y", space = "free_y")
ggsave(filename = paste0(path_to_inout, "ridge_kmer_cov_simple_", date, ".pdf"),
       plot = p4, device = "pdf", width = 10, height = 100, units = "cm")

p5 <- ggplot(data = df_denovo_non_zero, aes(x = log10(contig_kmer_cov), y = sample_id_kmer_cov_ord))
p5 <- p5 + geom_density_ridges(aes(fill = log10(total_num_of_contigs)), scale = 2, size = 0.5, stat = "binline", bins = 100)
p5 <- p5 + theme_ridges(font_size = 6, grid = FALSE)
p5 <- p5 + scale_fill_viridis_c()
p5 <- p5 + facet_grid(LibPool~., scales = "free_y", space = "free_y")
ggsave(filename = paste0(path_to_inout, "ridge_kmer_cov_LibPool_", date, ".pdf"),
       plot = p5, device = "pdf", width = 10, height = 100, units = "cm")

p6 <- ggplot(data = df_denovo_non_zero, aes(x = log10(contig_kmer_cov), y = sample_id_kmer_cov_ord))
p6 <- p6 + geom_density_ridges(aes(fill = log10(total_num_of_contigs)), scale = 2, size = 0.5, stat = "binline", bins = 100)
p6 <- p6 + theme_ridges(font_size = 6, grid = FALSE)
p6 <- p6 + scale_fill_viridis_c()
p6 <- p6 + facet_grid(genus_species~., scales = "free_y", space = "free_y")
ggsave(filename = paste0(path_to_inout, "ridge_kmer_cov_sp_", date, ".pdf"),
       plot = p6, device = "pdf", width = 10, height = 100, units = "cm")

df_denovo_nz_summary[which(df_denovo_nz_summary$total_num_of_contigs < 10),]



p61 <- ggplot(data = df_denovo_non_zero, aes(x = log10(contig_kmer_cov), y = sample_id_kmer_cov_ord))
p61 <- p61 + geom_density_ridges(aes(fill = log10(lensum_of_contigs)), scale = 2, size = 0.5, stat = "binline", bins = 100)
p61 <- p61 + theme_ridges(font_size = 6, grid = FALSE)
p61 <- p61 + scale_fill_viridis_c()
p61 <- p61 + facet_grid(LibPool~., scales = "free_y", space = "free_y")
ggsave(filename = paste0(path_to_inout, "ridge_kmer_cov_LibPool_lensumcolor_", date, ".pdf"),
       plot = p61, device = "pdf", width = 10, height = 100, units = "cm")



p7 <- ggplot()
p7 <- p7 + geom_histogram(data = df_denovo_nz_summary, mapping = aes(x = log10(total_num_of_contigs)), bins = 200)
p7 <- p7 + theme_minimal(base_size = 11)
#p7
ggsave(filename = paste0(path_to_inout, "hist_lg_contig_num_", date, ".pdf"),
       plot = p7, device = "pdf", width = 15, height = 10, units = "cm")

p8 <- ggplot()
p8 <- p8 + geom_histogram(data = df_denovo_nz_summary, mapping = aes(x = total_num_of_contigs), bins = 200)
p8 <- p8 + theme_minimal(base_size = 11)
#p8
ggsave(filename = paste0(path_to_inout, "hist_contig_num_", date, ".pdf"),
       plot = p8, device = "pdf", width = 15, height = 10, units = "cm")


#################################################

# Scatter av kmer cov-length vs av len
p21 <- ggplot(data = df_denovo_nz_summary,
             mapping = aes(x = log10(average_contig_len), y = log10(average_contig_kmer_cov)))
p21 <- p21 + geom_point(aes(color = LibPool, shape = NTC), size = 2)
p21 <- p21 + theme_minimal(base_size = 11)
p21

ggsave(filename = paste0(path_to_inout, "scatter_kmer_cov_vs_len_", date, ".pdf"),
       plot = p21, device = "pdf", width = 15, height = 15, units = "cm")

#### scatterplot matrix
cat(paste0(colnames(df_denovo_nz_summary), collapse = "', '"))


df_denovo_nz_summary$lg_total_num_of_contigs <- log10(df_denovo_nz_summary$total_num_of_contigs)
df_denovo_nz_summary$lg_lensum_of_contigs <- log10(df_denovo_nz_summary$lensum_of_contigs)
df_denovo_nz_summary$lg_average_contig_kmer_cov <- log10(df_denovo_nz_summary$average_contig_kmer_cov)
df_denovo_nz_summary$lg_num_qualtrim_kept_PE_reads <- log10(df_denovo_nz_summary$num_qualtrim_kept_PE_reads)
df_denovo_nz_summary$lg_total_num_of_PE_reads <- log10(df_denovo_nz_summary$total_num_of_PE_reads)

pscm1 <- ggpairs(df_denovo_nz_summary, columns = c('lg_total_num_of_PE_reads', 'lg_num_qualtrim_kept_PE_reads','prop_qualtrim_kept_PE_reads',
                                          'average_contig_len', 'lg_average_contig_kmer_cov',
                                          'lg_total_num_of_contigs', 'lg_lensum_of_contigs'), 
        ggplot2::aes(colour=LibPool)) 

pscm1 <- pscm1 + theme_minimal(base_size = 8)

ggsave(filename = paste0(path_to_inout, "scatter_matrix_", date, ".pdf"),
       plot = pscm1, device = "pdf", width = 40, height = 30, units = "cm")

###### only samples, that passed the trimming well with at least 1M reads


df_trim1 <- read_delim(file = "/full_path_to/wd/QC/analysis/Trimmomatic/trimmomatic_summary_CMVM_20211006.txt",
                      delim = "\t") 

OK_seq_and_trim <- df_trim1[which(df_trim1$flag_poor_seq_depth == 0 & df_trim1$flag_poor_read_qual == 0),]$sample_id

df_denovo_nz_summary_trimOK <- df_denovo_nz_summary[which(df_denovo_nz_summary$sample_id %in% OK_seq_and_trim),]

pscm11 <- ggpairs(df_denovo_nz_summary_trimOK, columns = c('lg_total_num_of_PE_reads', 'lg_num_qualtrim_kept_PE_reads','prop_qualtrim_kept_PE_reads',
                                                   'average_contig_len', 'lg_average_contig_kmer_cov',
                                                   'lg_total_num_of_contigs', 'lg_lensum_of_contigs'), 
                 ggplot2::aes(colour=LibPool)) 

pscm11 <- pscm11 + theme_minimal(base_size = 8)

ggsave(filename = paste0(path_to_inout, "scatter_matrix_seq_trimOK_", date, ".pdf"),
       plot = pscm11, device = "pdf", width = 50, height = 50, units = "cm")



## See correlation by batch

df_denovo_nz_summary_trimOK$post_extraction_Batch <- str_remove(as.character(df_denovo_nz_summary_trimOK$Batch), "\\.[[:digit:]]+")

df_denovo_nz_summary_trimOK[which(df_denovo_nz_summary_trimOK$post_extraction_Batch %in% c("13", "14", "15")),]$post_extraction_Batch <- "13_14_15"
df_denovo_nz_summary_trimOK[which(df_denovo_nz_summary_trimOK$post_extraction_Batch %in% c("16", "17", "18")),]$post_extraction_Batch <- "16_17_18"


paste0(unique(df_denovo_nz_summary_trimOK$post_extraction_Batch), collapse = "', '")

df_denovo_nz_summary_trimOK$post_extraction_Batch <- factor(df_denovo_nz_summary_trimOK$post_extraction_Batch, levels = c('1', '2', '3', 
                                                       '4', '5', '7',
                                                       '6', '8', '10', 
                                                       '11','12', 
                                                       '13_14_15',
                                                       '16_17_18'), 
                                                       ordered = T)

pscm2 <- ggpairs(df_denovo_nz_summary_trimOK, columns = c('lg_total_num_of_PE_reads', 'prop_qualtrim_kept_PE_reads',
                                                           'average_contig_len', 'lg_average_contig_kmer_cov',
                                                           'lg_total_num_of_contigs'), 
                  ggplot2::aes(colour=post_extraction_Batch)) 
pscm2 <- pscm2 + scale_color_ordinal()
pscm2 <- pscm2 + theme_minimal(base_size = 8)

ggsave(filename = paste0(path_to_inout, "scatter_matrix_seq_trimOK_by_post_extr_Batch_", date, ".pdf"),
       plot = pscm2, device = "pdf", width = 50, height = 50, units = "cm")



