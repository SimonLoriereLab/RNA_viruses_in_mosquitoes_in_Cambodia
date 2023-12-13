library(tidyverse)
library(seqinr)

N_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_N.txt", delim = "\t")
length(unique(N_df$contig_id))


G_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/contig_set_G.txt", delim = "\t") 

dmnd_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/lineage_complete.DIAMOND_hits_hmmscan_contigs.txt", delim = "\t",
                      guess_max = 17479)
                      
GN_contigs <- unique(c(G_df$contig_id, N_df$contig_id))


intersect(G_df$contig_id, N_df$contig_id)


GN_dmnd_df <- dmnd_df[which(dmnd_df$qseqid %in% GN_contigs),]


count = 0
for (q in unique(GN_dmnd_df$qseqid)){
  GN_dmnd_df_q <- GN_dmnd_df[which(GN_dmnd_df$qseqid == q),] %>% arrange(desc(bitscore))
  GN_dmnd_df_q_bb <- GN_dmnd_df_q[1,]
  if (count == 0){
    GN_dmnd_df_bb <- GN_dmnd_df_q_bb
  }else{
    GN_dmnd_df_bb <- bind_rows(GN_dmnd_df_bb, GN_dmnd_df_q_bb)
  }
  count = count + 1
}

intersect(G_df$contig_id, setdiff(GN_contigs, GN_dmnd_df_bb$qseqid))
intersect(N_df$contig_id, setdiff(GN_contigs, GN_dmnd_df_bb$qseqid))
N_contigs_wo_dmnd_hits <- setdiff(GN_contigs, GN_dmnd_df_bb$qseqid)
### 5 contigs without diamond hits are from the N set

### Select GN contigs for which best hits aren't viruses
GN_dmnd_df_bb_not_vir <- GN_dmnd_df_bb[which(!GN_dmnd_df_bb$superkingdom == "Viruses" | is.na(GN_dmnd_df_bb$superkingdom)),]
## make sure they are from N
intersect(GN_dmnd_df_bb_not_vir$qseqid, G_df$contig_id)
intersect(GN_dmnd_df_bb_not_vir$qseqid, N_df$contig_id)


N_hmm_df<- read_delim("/full_path_to/wd/RdRp_scan/analysis/N_df.txt", delim = "\t")

GN_dmnd_df_bb_not_vir1 <- left_join(GN_dmnd_df_bb_not_vir, N_hmm_df, by = c("qseqid"="contig_id"))


colnames(GN_dmnd_df_bb_not_vir1)


p1 <- ggplot()
p1 <- p1 + geom_point(data = GN_dmnd_df_bb_not_vir1, 
                      mapping = aes(x = bitscore, y = score_full_seq))


plot(p1)

##### Extract all GN ORFs and rerun blastp


GN_dmnd_df_bb_not_vir1$rdrp_scan_orf_id

all_hmm_df <- read_delim("/full_path_to/wd/RdRp_scan/2022_07_21_rdrp_results/code_1/hmmscan_tblout_best_AB.txt", delim = "\t")
all_hmm_df <- all_hmm_df %>% 
  mutate_all(.funs = str_replace_all, ",", ".") %>%
  mutate_at(c("evalue_full_seq", "score_full_seq", "bias_full_seq", 
              "evalue_best_dom", "score_best_dom", "bias_best_dom"),
            .funs = as.numeric)

all_hmm_df$contig_id <- str_remove(all_hmm_df$rdrp_scan_orf_id, "_[[:digit:]]+$")

GN_hmm_df <- all_hmm_df[which(all_hmm_df$contig_id %in% GN_contigs),]

GN_hmm_df$contig_set <- "N"
GN_hmm_df[which(GN_hmm_df$contig_id %in% G_df$contig_id),]$contig_set <- "G"

unique(GN_hmm_df$rdrp_scan_orf_id)

orf_nr_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/all_emboss_getorfs_nr.txt", delim = "\t")

GN_hmm_df1 <- left_join(GN_hmm_df, orf_nr_df, by = c("contig_id", "rdrp_scan_orf_id" = "emboss_getorf_id"))
GN_hmm_df1$orf_len <- abs(as.numeric(GN_hmm_df1$orf_end) - as.numeric(GN_hmm_df1$orf_start))
GN_hmm_df1_300aa <- GN_hmm_df1[which(GN_hmm_df1$orf_len >= 900),]

length(unique(GN_hmm_df1_300aa$rdrp_scan_orf_id))
length(unique(GN_hmm_df1_300aa$contig_id))




write.table(x = GN_hmm_df1_300aa, file = "/full_path_to/wd/RdRp_scan/analysis/GN_hmmscan_df.txt", append = F,quote = F,sep = "\t",
            row.names = F, col.names = T)


# orfs_nr <- read.fasta("/full_path_to/wd/RdRp_scan/analysis/all_samples_renamed_contigs_redundancy_removed_code1_aa_nr.fasta", strip.desc = T)
# 
# orfs_nr_rdrp <- orfs_nr[GN_hmm_df1_300aa$rdrp_scan_orf_id]

# write.fasta(sequences = orfs_nr_rdrp, 
#             names = names(orfs_nr_rdrp), 
#             file.out = "/full_path_to/wd/RdRp_scan/analysis/all_samples_renamed_contigs_redundancy_removed_code1_aa_nr_rdrp.fasta")
