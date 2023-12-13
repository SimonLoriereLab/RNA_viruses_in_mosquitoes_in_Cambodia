#### Run on MAESTRO HPC

library(tidyverse)
library(lubridate)

#### Find all paths
path_to_sterr <- "/full_path_to/wd/auto_metaviromics_jobserrors_trim/"
# path_to_sterr <- "/full_path_to/wd/auto_metaviromics_jobserrors/"

#### Output path
path_to_output <- "/full_path_to/wd/QC/analysis/Trimmomatic/"

#### Function to extract info from single error file
ext_readnum_fun <- function(path_to_err_file){
  print(path_to_err_file)
  err_lines <- read_lines(file = path_to_err_file, num_threads = 95)
  mtime <- ymd_hms(file.info(path_to_err_file)$mtime, tz = "CET")
  err_fname <- rownames(file.info(path_to_err_file))
  for (l in err_lines){
    if(str_detect(l, "^ -phred33")){
      sample_id <- str_extract(l, "(?<=/)[[:alnum:]]+(?=_R1.fastq.gz?)")
    }else if (str_detect(l, "^Input Read Pairs")){
      input_read_pairs <- as.numeric(str_remove_all(str_extract(l, "(?<=^Input Read Pairs:)[\\s]+[[:digit:]]+"), "[\\s]*"))
      both_surv <- as.numeric(str_remove_all(str_extract(l, "(?<=Both Surviving:)[\\s]+[[:digit:]]+"), "[\\s]*"))
      for_only_surv <- as.numeric(str_remove_all(str_extract(l, "(?<=Forward Only Surviving:)[\\s]+[[:digit:]]+"), "[\\s]*"))
      rev_only_surv <- as.numeric(str_remove_all(str_extract(l, "(?<=Reverse Only Surviving:)[\\s]+[[:digit:]]+"), "[\\s]*"))
      dropped <- as.numeric(str_remove_all(str_extract(l, "(?<=Dropped:)[\\s]+[[:digit:]]+"), "[\\s]*"))
      # check
      both_surv + for_only_surv + rev_only_surv + dropped == input_read_pairs
      total_reads <- input_read_pairs*2
      break
    }
  }
  if (exists("sample_id")){
    return(tibble(sample_id, mtime, total_reads, input_read_pairs, both_surv, for_only_surv, rev_only_surv, dropped, err_fname))
  }else{
    print("no trimming in this error file")
    return(NULL)
  }
}

#### Run scraping
count = 0
for (i in path_to_sterr){
  #print(i)
  for (j in list.files(path = i, pattern = ".err$", all.files = FALSE, full.names = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)){
    #print(j)
    path_to_err_file <- paste0(i, j)
    print(path_to_err_file)
    read_num_single_sample_df <- ext_readnum_fun(path_to_err_file)
    if (is.null(read_num_single_sample_df)){
      next
    }
    count = count + 1
    if (count == 1){
      read_num_all_df <- read_num_single_sample_df
    }else{
      read_num_all_df <- bind_rows(read_num_all_df,read_num_single_sample_df)
    }
  }
}


read_num_all_df$total_num_of_PE_reads <- read_num_all_df$input_read_pairs*2

read_num_all_df$num_qualtrim_kept_PE_reads <- read_num_all_df$both_surv*2

read_num_all_df$prop_qualtrim_kept_PE_reads <- read_num_all_df$num_qualtrim_kept_PE_reads/read_num_all_df$total_num_of_PE_reads

read_num_all_df1 <- read_num_all_df %>% select(sample_id, err_fmtime = mtime, err_fname, 
                                               total_num_of_PE_reads, num_qualtrim_kept_PE_reads, prop_qualtrim_kept_PE_reads, for_only_surv, rev_only_surv, dropped)

date <- str_remove_all(Sys.Date(), "-")

# write.table(x = read_num_all_df1, file = paste0(path_to_output, "trimmomatic_scraped_stats_CMVM_", date, ".txt"),
#             append = F, quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = read_num_all_df1, file = paste0(path_to_output, "trimmomatic_scraped_stats_CMVM_", date, "_1.txt"),
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


