library(tidyverse)

### Correct the directories if needed
main_directory <- "/full_path_to/wd/rRNA_depletion_verification" ### from MAESTRO HPC

# sample_info <- read.delim(paste0(main_directory,
#                                  "/metadata/CMVM - priority annotated full sample list 26072021.tsv"))
input_directory <- paste0(main_directory,
                          "/analysis/Aedes_Culex_depleted/mapping_rRNA")
output_directory <- paste0(main_directory,
                            "/analysis/Aedes_Culex_depleted/mapping_analysis_rRNA")
dir.create(output_directory)


files_summary = list.files(input_directory, pattern="*mapsum.txt.gz")
#print(files_summary)
count_sum = 1
for (n in 1:length(files_summary)) {
  #import each dataframe according to the list set above
  #print(files_summary[n])
  file_inf <- file.info(paste0(input_directory,"/", files_summary[n]))
  if (file_inf$size != 0) {
    print(count_sum)
    #cat("\nnot empty file\n")
    df_sum_temp <- read.delim(paste0(input_directory,"/", files_summary[n]))
    if (count_sum == 1){
      df_sum_final <- df_sum_temp
      count_sum = count_sum + 1
    }else{
      df_sum_final <- rbind(df_sum_final, df_sum_temp)
      count_sum = count_sum + 1
    }
  }
}

write.table(x = df_sum_final, file = gzfile(paste0(output_directory,"/all_multimap_mapsum_rRNA.txt.gz")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

