library(tidyverse)
library(seqinr)
library(rentrez)

dir.create("/full_path_to/wd/RdRp_scan/analysis/too_short_contigs_analysis/")

L_set_df2_be <- read_delim("/full_path_to/wd/RdRp_scan/analysis/L_df_hmmscan.LCA.BestEvalue_hit.txt", 
                         delim = "\t")
### Not necessary to consider G for this
L_set_df2_be_N <- L_set_df2_be[which(L_set_df2_be$contig_set == "N" & 
                                       L_set_df2_be$superkingdom == "Viruses"),]


### Species LCA available
L_set_df2_be_N_species_LCA <- L_set_df2_be_N[which(!is.na(L_set_df2_be_N$species)),]
L_set_df2_be_N_species_LCA$taxon_rank_genome_segment_size <- "species"
species_LCA_taxids <- sort(unique(L_set_df2_be_N_species_LCA$tax_id))
### Genus LCA available
L_set_df2_be_N_genus_LCA <- L_set_df2_be_N[which(!is.na(L_set_df2_be_N$genus) & 
                                                !L_set_df2_be_N$contig_id %in% unique(L_set_df2_be_N_species_LCA$contig_id)),]
L_set_df2_be_N_genus_LCA$taxon_rank_genome_segment_size <- "genus"
genus_LCA_taxids <- sort(unique(L_set_df2_be_N_genus_LCA$tax_id))

### Family LCA available
L_set_df2_be_N_family_LCA <- L_set_df2_be_N[which(!is.na(L_set_df2_be_N$family) & 
                                                   !L_set_df2_be_N$contig_id %in% unique(L_set_df2_be_N_species_LCA$contig_id) &
                                                    !L_set_df2_be_N$contig_id %in% unique(L_set_df2_be_N_genus_LCA$contig_id)),]

L_set_df2_be_N_family_LCA$taxon_rank_genome_segment_size <- "family"
family_LCA_taxids <- sort(unique(L_set_df2_be_N_family_LCA$tax_id))


LuN_LCA_genus_family <- bind_rows(L_set_df2_be_N_genus_LCA, L_set_df2_be_N_family_LCA)



LuN_LCA_shallow <- bind_rows(L_set_df2_be_N_species_LCA, L_set_df2_be_N_genus_LCA, L_set_df2_be_N_family_LCA)


shallow_LCA_taxids <- LuN_LCA_shallow %>% 
  select(tax_id, taxon_rank_genome_segment_size) %>% distinct()

LuN_LCA_deep <- L_set_df2_be_N[which(!L_set_df2_be_N$contig_id %in% LuN_LCA_shallow$contig_id),]

LuN_LCA_deep$taxon_rank_genome_segment_size <- "species_be"

deep_LCA_taxids <- LuN_LCA_deep %>% 
  select(tax_id = tax_id_be, taxon_rank_genome_segment_size) %>% distinct()

LuN_LCA_shallow_and_deep <- bind_rows(LuN_LCA_shallow, LuN_LCA_deep)

LCA_taxids_shallow_and_deep <- bind_rows(shallow_LCA_taxids, deep_LCA_taxids) %>% arrange(tax_id)
cat(paste0(LCA_taxids_shallow_and_deep$tax_id, collapse = ","))

### Annotate expected genome/segment sizes for some well-studied taxons
write.table(x = LCA_taxids_shallow_and_deep[which(LCA_taxids_shallow_and_deep$taxon_rank_genome_segment_size %in% c("family", "genus")),],
            file = "/full_path_to/wd/RdRp_scan/analysis/too_short_contigs_analysis/family_genus_for_genome_segment_size_annotation.txt", 
            sep = "\t", append = F, quote = F, row.names = F, col.names = T)

##### Extract length and titles for all sequences in these taxons
extract_accessions_and_slen_by_taxid <- function(input_tax_ids_vector, output_directory, ncbi_api_key){
  
  dir.create(output_directory)
  
  count = 0
  for (itxid in input_tax_ids_vector){
    print("Input taxid:")
    print(itxid)
    for (i in seq(from = 100, to = 1000000, by = 100)){
      print (i)
      tsearch <- entrez_search(db="taxonomy", term= paste0("txid", itxid,"[Subtree]"), retmax = i,  api_key = ncbi_api_key)
      tsearch$ids
      if (i > length(tsearch$ids)){
        break
      }
    }
    
    summary_partition_size <- 100
    n_summary_partitions <- (length(tsearch$ids)%/%summary_partition_size + 1)
    for (summary_partition in (1:n_summary_partitions)){
      
        print(paste0("processing summary_partition ", summary_partition, "/", n_summary_partitions))
        
        if (summary_partition == 1){
          txid_min <- 1
        }else{
          txid_min <- txid_max + 1
        }
        
        if (summary_partition == n_summary_partitions){
          txid_max <- length(tsearch$ids)
        }else{
          txid_max <- summary_partition*summary_partition_size
        }
        
        tsum <- entrez_summary(db="taxonomy", id=tsearch$ids[txid_min:txid_max], 
                               api_key = ncbi_api_key)

        
        if(length(tsearch$ids[txid_min:txid_max]) > 1){
          tax_df <- data.frame(do.call(rbind,tsum))
        }else{
          tax_df <- data.frame(tsum)
        }

        tax_df1 <- tax_df[which(tax_df$status == "active" & 
                                  tax_df$rank %in% c("species", "isolate") &
                                  tax_df$division == "viruses" &
                                  tax_df$genbankdivision == "Viruses"),]
        tax_ids <- as.vector(unlist(tax_df1$taxid))
        
        
        df_nt <- NULL
        nt_count = 0
        if (length(tax_ids)>0){
          for (taxid in tax_ids){
            #taxid <- "11292"
            cat("\n\n")
            print(paste0("Processing tax ID ", taxid))
            n_unprocessed_taxids <- (length(tax_ids)-match(taxid, tax_ids))
            print(paste0(n_unprocessed_taxids, " tax IDs remain to be processed"))
            
            ntsearch0 <- entrez_search(db="nuccore", term = paste0("txid", taxid, "[taxid]"), 
                                       api_key = ncbi_api_key, retmax=0)
            total_n_ids <- ntsearch0$count
            print("total seq IDs to process")
            print(total_n_ids)
            ntsearch <- entrez_search(db="nuccore", term = paste0("txid", taxid, "[taxid]"), 
                                      api_key = ncbi_api_key, retmax=total_n_ids)
            cat("\n\n")
            
            if (length(ntsearch$ids)>0){ 
              #### partitioning processing of IDs, because entrez_summary cannot take large numbers of IDs at once
              partition_size <- 100
              n_partitions <- (ntsearch0$count%/%partition_size + 1)
              
              if (n_partitions > 3){
                
                ### pick first, last and 2 random in the middle
                partitions_to_analyze <- c(1, sample(2:(n_partitions-1), 2), n_partitions)
              }else{
                partitions_to_analyze <- c(1:n_partitions)
              }
              
              
              for (partition in partitions_to_analyze){
                
                print(paste0("processing partition ", partition, "/", length(partitions_to_analyze)))
                
                if (partition == 1){
                  id_min <- 1
                }else{
                  id_min <- (partition-1)*partition_size + 1
                }
                if (partition == n_partitions){
                  id_max <- total_n_ids
                }else{
                  id_max <- partition*partition_size
                }
                nt_sum <- entrez_summary(db="nuccore", id = ntsearch$ids[id_min:id_max], 
                                         retmode = "json", always_return_list = T, 
                                         api_key = ncbi_api_key, retmax=1000)
                
                #### Process each partition of IDs. That is useful not to call each id summary separately.
                for(j in 1:length(nt_sum)){
                  #### j = 2 #########################################
                  ntdfgi <- as_tibble(nt_sum[[j]][1:29]) %>% select(nt_uid=uid, nt_accession = caption, nt_organism = organism, nt_title = title, nt_slen =slen)
                  ntdfgi$tax_id <- taxid
                  ntdfgi$total_number_of_sequences <- total_n_ids
                  if (nt_count == 0){
                    df_nt <- ntdfgi
                  }else{
                    df_nt <- bind_rows(df_nt, ntdfgi)
                  }
                  nt_count = nt_count +1
                }
              }
            }else{
              print("nt entry not found by ncbi txid")
            }
          }
      }
      if (!is.null(df_nt)){
        df_nt$input_tax_id <- itxid
        if (count == 0){
          df_nt_all <- df_nt
        }else{
          df_nt_all <- bind_rows(df_nt_all, df_nt)
        }
        count = count + 1
      }
    }
  }
  
  return(df_nt_all)

}



output_directory <- "/full_path_to/wd/RdRp_scan/analysis/too_short_contigs_analysis/NCBI_API_extract_accessions_and_slen_by_taxid"
ncbi_api_key <- "b846efefdc1e057fae015142125e6842d108"
LCA_taxids_shallow_and_deep$tax_id


# taxid_lengths_titles_df <- extract_accessions_and_slen_by_taxid(
#   input_tax_ids_vector=LCA_taxids_shallow_and_deep[which(LCA_taxids_shallow_and_deep$taxon_rank_genome_segment_size %in%
#                                                            c("species", "species_be")),]$tax_id,
#                                      output_directory=output_directory,
#                                      ncbi_api_key=ncbi_api_key)
# 
# 
# ##### Verify if there's any mismatch of tax_id and input tax_id
# sum(taxid_lengths_titles_df$tax_id == taxid_lengths_titles_df$input_tax_id) == nrow(taxid_lengths_titles_df)
# 
# write.table(x=taxid_lengths_titles_df, file = paste0(output_directory, "/species_NCBI_titles_and_lengths.txt"),
#             col.names = T, row.names = F, append = F, quote = F, sep = "\t")


## Keep only complete entries
taxid_lengths_titles_df <- read_delim( paste0(output_directory, "/species_NCBI_titles_and_lengths.txt"), delim = "\t")
taxid_lengths_titles_df1 <- taxid_lengths_titles_df[which(str_detect(taxid_lengths_titles_df$nt_title, ", complete")),]

## "2585031" this is Riboviria sp., which is in fact not a single species, so I'll ignore it.

taxid_lengths_titles_df2 <- taxid_lengths_titles_df1[which(!taxid_lengths_titles_df1$input_tax_id == "2585031"),]

unique(taxid_lengths_titles_df$input_tax_id)
unique(taxid_lengths_titles_df1$input_tax_id)
unique(taxid_lengths_titles_df2$input_tax_id)

taxid_lengths_titles_df3 <- taxid_lengths_titles_df2 %>% 
  select(input_tax_id, nt_slen, nt_title, nt_organism, nt_accession) %>% distinct()
taxid_lengths_titles_df3$tax_id_nt_organism <- paste0(taxid_lengths_titles_df3$input_tax_id, " - ",taxid_lengths_titles_df3$nt_organism)
LuN_LCA_shallow_and_deep$taxon_rank_genome_segment_size
LuN_LCA_shallow_and_deep$tax_id




LuN_LCA_shallow_and_deep$length_analysis_input_tax_id <- LuN_LCA_shallow_and_deep$tax_id

LuN_LCA_shallow_and_deep[which(
  LuN_LCA_shallow_and_deep$taxon_rank_genome_segment_size == "species_be"), 
  ]$length_analysis_input_tax_id <- LuN_LCA_shallow_and_deep[
    which(LuN_LCA_shallow_and_deep$taxon_rank_genome_segment_size == "species_be"), 
    ]$tax_id_be

taxid_lengths_titles_df3_LuN <- left_join(taxid_lengths_titles_df3, 
                                          LuN_LCA_shallow_and_deep,
                                          by = c("input_tax_id"="length_analysis_input_tax_id")) %>% 
  distinct()


taxid_lengths_titles_df3_LuN$contig_len <- as.numeric(
  str_extract(taxid_lengths_titles_df3_LuN$contig_id,
  "(?<=length\\_)[[:digit:]]+"))



contig_lens <- taxid_lengths_titles_df3_LuN %>% select(input_tax_id, contig_id, contig_len, nt_organism) %>% 
  distinct()

contig_lens$tax_id_nt_organism <- paste0(contig_lens$input_tax_id, " - ",contig_lens$nt_organism)


##### Plot size distributions and actual sizes

taxid_lengths_titles_df4 <- full_join(taxid_lengths_titles_df3, contig_lens) 
count2 = 0
for (txo in unique(taxid_lengths_titles_df4$tax_id_nt_organism)){
  tdf <- taxid_lengths_titles_df4[which(taxid_lengths_titles_df4$tax_id_nt_organism ==  txo),] %>% 
    arrange(contig_len)
  tdf1 <- tdf[which(tdf$nt_slen == min(tdf$nt_slen)),][1,]
  tdf1$diff <- tdf1$nt_slen - tdf1$contig_len
  tdf1$perc <- floor((tdf1$contig_len/tdf1$nt_slen)*100)
  tdf2 <- tdf1 %>% select(tax_id_nt_organism, input_tax_id, diff, perc, contig_id, contig_len)
  if (count2 == 0){
    ordering_df <- tdf2
  }else{
    ordering_df <- bind_rows(ordering_df, tdf2)
  }
  count2 = count2 + 1
}

ordering_df1 <- ordering_df %>% arrange(desc(diff))



taxid_lengths_titles_df3$tax_id_nt_organism <- factor(taxid_lengths_titles_df3$tax_id_nt_organism, 
                                                      levels = unique(ordering_df1$tax_id_nt_organism), 
                                                      ordered = T)
contig_lens$tax_id_nt_organism <- factor(contig_lens$tax_id_nt_organism, 
                                                      levels = unique(ordering_df1$tax_id_nt_organism), 
                                         ordered = T)
ordering_df1$tax_id_nt_organism <- factor(ordering_df1$tax_id_nt_organism, 
                                         levels = unique(ordering_df1$tax_id_nt_organism), 
                                         ordered = T)

p1 <- ggplot()
p1 <- p1 + geom_jitter(data = taxid_lengths_titles_df3, 
                      mapping = aes(y = tax_id_nt_organism, x = nt_slen), 
                      shape = 19, color = "gray50")

p1 <- p1 + geom_point(data = contig_lens, 
                       mapping = aes(y = tax_id_nt_organism, x = contig_len), 
                       shape = 0, color  = "red", size = 4)
p1 <- p1 + geom_text(data = ordering_df1, 
                      mapping = aes(y = tax_id_nt_organism, x = contig_len, label = perc),
                     color  = "black", size =2)
p1 <- p1 + theme_minimal(base_size = 9)
p1 <- p1 + theme(axis.text.y = element_text(colour = ifelse(test = ordering_df1$input_tax_id %in% deep_LCA_taxids$tax_id, 
                                                            yes = "blue", no = "black")))


print(p1)

ggsave(plot = p1, 
       filename = paste0(output_directory, "/contig_size_comparison_to_NT_available_references_temp_____.pdf"),
       device = "pdf",width = 20, height = 29, dpi = 300, units = "cm")



###### Contigs to remove due to very small size comparing to whats expected in the database
###### Check similarity level before removal



contigs_to_remove_df <- L_set_df2_be_N[which(L_set_df2_be_N$tax_id %in% 
                                               c("1922677", "2662141", "1922360",
                                                                          "936308", "2847850", "1923779",
                                                                          "2010280", "1923775")),] %>% 
  select(tax_id, contig_id, pident_be, length_be, stitle_be, tax_id_be)



contigs_to_remove_df_be <- L_set_df2_be_N[which(L_set_df2_be_N$tax_id_be %in% 
                                               c("1923161", "2749925", "1923078", "2899248", "2899247", "2587545") & 
                                                 !L_set_df2_be_N$contig_id %in% c("F1628_NODE_3_length_6523_cov_75.503556",
                                                                                 "F1015_NODE_34_length_4438_cov_44.128451",
                                                                                 "F1194_NODE_83_length_4344_cov_522.683143")),] %>% 
  select(tax_id, contig_id, pident_be, length_be, stitle_be, tax_id_be)


contigs_removed_too_short_df <- bind_rows(contigs_to_remove_df, contigs_to_remove_df_be)

contig_ids_to_remove <- unique(contigs_removed_too_short_df$contig_id)

contig_lens$too_short <- "no" 
contig_lens[which(contig_lens$contig_id %in% contig_ids_to_remove),]$too_short <- "yes"


write.table(x = contig_ids_to_remove,
            file = "/full_path_to/wd/RdRp_scan/analysis/too_short_contigs_analysis/too_short_contigs_to_remove_from_L_species_level_LCA_or_BestEvalue_hit_above_family_level.txt", 
            sep = "\t", append = F, quote = F, row.names = F, col.names = F)





p2 <- ggplot()
p2 <- p2 + geom_jitter(data = taxid_lengths_titles_df3, 
                       mapping = aes(y = tax_id_nt_organism, x = nt_slen), 
                       shape = 19, color = "gray50")

p2 <- p2 + geom_point(data = contig_lens, 
                      mapping = aes(y = tax_id_nt_organism, x = contig_len, shape = too_short), 
                      color  = "red", size = 4)
p2 <- p2 + scale_shape_manual(breaks = c("yes" , "no"), values = c(4, 0))
p2 <- p2 + geom_text(data = ordering_df1, 
                     mapping = aes(y = tax_id_nt_organism, x = contig_len, label = perc),
                     color  = "black", size = 2)
p2 <- p2 + theme_minimal(base_size = 9)
p2 <- p2 + theme(axis.text.y = element_text(colour = ifelse(test = ordering_df1$input_tax_id %in% deep_LCA_taxids$tax_id, 
                                                            yes = "blue", no = "black")))



print(p2)

ggsave(plot = p2, 
       filename = paste0(output_directory, "/contig_removal_by_size_species_level.pdf"),
       device = "pdf",width = 20, height = 29, dpi = 300, units = "cm")


check_df <- taxid_lengths_titles_df3[which(taxid_lengths_titles_df3$input_tax_id == 2587545),]

check_df1 <- contig_lens[which(contig_lens$input_tax_id == 2587545),]
check_df2 <- ordering_df1[which(ordering_df1$input_tax_id == 2587545),]
########## 



###########################################################################################
###########################################################################################
### Annotate expected genome/segment sizes for some well-studied taxons
family_genus_genome_size_df <- read_delim("/full_path_to/wd/RdRp_scan/analysis/too_short_contigs_analysis/CMVM_bioinformatics - genome_segment_size.tsv", 
           delim = "\t")

LuN_LCA_genus_family1 <- left_join(LuN_LCA_genus_family, family_genus_genome_size_df, 
                                   by = c("tax_id", "taxon_rank_genome_segment_size"))

LuN_LCA_genus_family1$size

LuN_LCA_genus_family1$contig_len <- as.numeric(
  str_extract(LuN_LCA_genus_family1$contig_id,
              "(?<=length\\_)[[:digit:]]+"))

LuN_LCA_genus_family1$LCA <- LuN_LCA_genus_family1$family
LuN_LCA_genus_family1[which(LuN_LCA_genus_family1$taxon_rank_genome_segment_size == "genus"), ]$LCA <- LuN_LCA_genus_family1[which(LuN_LCA_genus_family1$taxon_rank_genome_segment_size == "genus"), ]$genus
LuN_LCA_genus_family1$tax_id_LCA <- paste0(LuN_LCA_genus_family1$tax_id, " - ", LuN_LCA_genus_family1$LCA)

contig_lens1 <- LuN_LCA_genus_family1 %>% select(tax_id, contig_id, contig_len, LCA, tax_id_LCA) %>% 
  distinct()

##### Plot size distributions and actual sizes

LuN_LCA_genus_family2 <- full_join(LuN_LCA_genus_family1, contig_lens1) 
LuN_LCA_genus_family2$tax_id_LCA

count3 = 0
for (txo in unique(LuN_LCA_genus_family2$tax_id_LCA)){
  tdf <- LuN_LCA_genus_family2[which(LuN_LCA_genus_family2$tax_id_LCA ==  txo),] %>% 
    arrange(contig_len)
  tdf1 <- tdf[which(tdf$size == min(tdf$size)),][1,]
  tdf1$diff <- tdf1$size - tdf1$contig_len
  tdf1$perc <- floor((tdf1$contig_len/tdf1$size)*100)
  tdf2 <- tdf1 %>% select(tax_id_LCA, tax_id, diff, perc, contig_id, contig_len)
  if (count3 == 0){
    family_genus_ordering_df <- tdf2
  }else{
    family_genus_ordering_df <- bind_rows(family_genus_ordering_df, tdf2)
  }
  count3 = count3 + 1
}

family_genus_ordering_df1 <- family_genus_ordering_df %>% arrange(desc(diff))



LuN_LCA_genus_family1$tax_id_LCA <- factor(LuN_LCA_genus_family1$tax_id_LCA, 
                                                      levels = unique(family_genus_ordering_df1$tax_id_LCA), 
                                                      ordered = T)
contig_lens1$tax_id_LCA <- factor(contig_lens1$tax_id_LCA, 
                                         levels = unique(family_genus_ordering_df1$tax_id_LCA), 
                                         ordered = T)
family_genus_ordering_df1$tax_id_LCA <- factor(family_genus_ordering_df1$tax_id_LCA, 
                                          levels = unique(family_genus_ordering_df1$tax_id_LCA), 
                                          ordered = T)

p3 <- ggplot()
p3 <- p3 + geom_jitter(data = LuN_LCA_genus_family1, 
                       mapping = aes(y = tax_id_LCA, x = size), 
                       shape = 19, color = "gray50")

p3 <- p3 + geom_point(data = contig_lens1, 
                      mapping = aes(y = tax_id_LCA, x = contig_len), 
                      shape = 0, color  = "red", size = 4)
p3 <- p3 + geom_text(data = family_genus_ordering_df1, 
                     mapping = aes(y = tax_id_LCA, x = contig_len, label = perc),
                     color  = "black", size =2)
p3 <- p3 + theme_minimal(base_size = 9)


print(p3)

ggsave(plot = p3, 
       filename = "/full_path_to/wd/RdRp_scan/analysis/too_short_contigs_analysis/contig_size_comparison_to_family_genus_known_size_temp_____.pdf",
       device = "pdf",width = 20, height = 29, dpi = 300, units = "cm")


contigs_family_genus_too_short <- c("F1611_NODE_535_length_1459_cov_2221.661681", "F1216_NODE_2_length_2921_cov_420.068388", "F1260_NODE_533_length_1135_cov_3.767593")

contig_lens1$too_short <- "no" 
contig_lens1[which(contig_lens1$contig_id %in% contigs_family_genus_too_short),]$too_short <- "yes"




p4 <- ggplot()
p4 <- p4 + geom_jitter(data = LuN_LCA_genus_family1, 
                       mapping = aes(y = tax_id_LCA, x = size), 
                       shape = 19, color = "gray50")

p4 <- p4 + geom_point(data = contig_lens1, 
                      mapping = aes(y = tax_id_LCA, x = contig_len, shape = too_short), 
                      color  = "red", size = 4)
p4 <- p4 + scale_shape_manual(breaks = c("yes" , "no"), values = c(4, 0))
p4 <- p4 + geom_text(data = family_genus_ordering_df1, 
                     mapping = aes(y = tax_id_LCA, x = contig_len, label = perc),
                     color  = "black", size = 2)
p4 <- p4 + theme_minimal(base_size = 9)



print(p4)

ggsave(plot = p4, 
       filename = "/full_path_to/wd/RdRp_scan/analysis/too_short_contigs_analysis/contig_removal_by_size_family_genus_level_____.pdf",
       device = "pdf",width = 20, height = 29, dpi = 300, units = "cm")





###### Remove all the contigs that were deemed too short from all taxonomic levels (all analyses above)
### Make set R: not too short L ready for first phylogeny
R_df <- L_set_df2_be[which(!L_set_df2_be$contig_id %in% c(contig_ids_to_remove,
                                                          contigs_family_genus_too_short)),]
write.table(x = R_df,
            file = "/full_path_to/wd/RdRp_scan/analysis/R_df_hmmscan.LCA.BestEvalue_hit_____.txt", 
            sep = "\t", append = F, quote = F, row.names = F, col.names = T)





