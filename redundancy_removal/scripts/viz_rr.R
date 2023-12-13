library(tidyverse)


# headers <- read_delim(file = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/headers_all_samples_renamed_contigs.fasta", 
#                                   num_threads = 16, delim = "\t", col_names = F)
# 
# contig_len <- as.numeric(str_extract(headers$X1, "(?<=length_)[[:digit:]]+"))
# 
# full_contig_df <- tibble(contig_names = headers$X1, contig_len)
# 
# write.table(x = full_contig_df, 
#             file = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/all_contig_names_len.txt", 
#             append = F,quote = F, sep = "\t",row.names = F,col.names = T)

full_contig_df <- read_delim(file = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/all_contig_names_len.txt",
                      num_threads = 16, delim = "\t")


nrow(full_contig_df)/10^6

# nr_contigs <- read_delim(file = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/headers_all_samples_renamed_contigs_redundancy_removed.fasta", 
#                                     num_threads = 16, delim = "\t", col_names = F)
# 
# 
# nr_contig_df <- full_contig_df[which(full_contig_df$contig_names %in% nr_contigs$X1),]
# 
# write.table(x = nr_contig_df, 
#             file = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/nr_contig_names_len.txt", 
#             append = F,quote = F, sep = "\t",row.names = F,col.names = T)

nr_contig_df <- read_delim(file = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/nr_contig_names_len.txt",
                      num_threads = 16, delim = "\t")

nrow(nr_contig_df)/10^6

# p0 <- ggplot()
# p0 <- p0 + geom_histogram(data = full_contig_df, aes(x = log10(contig_len)), bins = 100, color = "gray")
# 
# ggsave(plot = p0, 
#        filename = "/Volumes/GEVA/users/Artem/CMVM_lp0_20_data/redundancy_removal/analysis/extracted_renamed_contigs/hist_all_contig_len.pdf",
#        device = "pdf", width = 10, height = 5, units = "cm")

# p1 <- ggplot()
# p1 <- p1 + geom_density(data = full_contig_df, aes(x = log10(contig_len)), color = "grey",fill = "grey", alpha = 0.5)
# p1 <- p1 + geom_density(data = nr_contig_df, aes(x = log10(contig_len)), color = "black",fill = "black", alpha = 0.5)
# p1 <- p1 + theme_minimal(base_size = 14, base_family = "Helvetica")

# ggsave(plot = p1, 
#        filename = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/density_contig_len.pdf",
#        device = "pdf", width = 20, height = 10, units = "cm")

# p1 <- ggplot()
# p1 <- p1 + geom_density(data = full_contig_df, aes(x = log10(contig_len)), color = "grey",fill = "grey", alpha = 0.1)
# p1 <- p1 + geom_density(data = nr_contig_df, aes(x = log10(contig_len)), color = "black",fill = "black", alpha = 0.1)
# p1 <- p1 + theme_minimal(base_size = 14, base_family = "Helvetica")

# ggsave(plot = p1, 
#        filename = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/density_contig_len1.pdf",
#        device = "pdf", width = 20, height = 10, units = "cm")

p1 <- ggplot()
p1 <- p1 + geom_density(data = full_contig_df, aes(x = log10(contig_len)), color = "grey",fill = "grey", alpha = 1)
p1 <- p1 + geom_density(data = nr_contig_df, aes(x = log10(contig_len)), color = "firebrick3",fill = "firebrick3", alpha = 0.1)
p1 <- p1 + xlab("Contig length, log10nt")
p1 <- p1 + theme_classic(base_size = 8)

ggsave(plot = p1, 
       filename = "/full_path_to/wd/redundancy_removal/analysis/extracted_renamed_contigs/density_contig_len_v2.pdf",
       device = "pdf", width = 7, height = 4, units = "cm")

