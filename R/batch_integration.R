# batch integration
library(tidyverse)
library(sva)
library(CRISPRcleanRatSquared)
library(reshape2)
setwd("/group/iorio/lucia/datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/")
source("batch_integration_functions.R")

lib_name <- c("COLO1", "COLO2", "COLO3", 
              "BRCA1", "BRCA2", "BRCA3")
data <- library <- list()

for (i in 1:length(lib_name)) {
  
  curr_lib <- lib_name[i]
  print(curr_lib)
  data[[i]] <- read_csv(sprintf("%s_FINAL_EXACT_logFC_sgRNA_scaled.txt", curr_lib))  
  
  # load lib
  if (grepl("COLO", curr_lib)) {
    part1 <- "COREAD"
  }else{
    part1 <- "BRCA"
  }
  part2 <- stringr::str_sub(curr_lib, start = 5, end = 5)
  lib_file_name <- sprintf("LIBS/ENCORE_GI_%s_Library_%s.txt", part1, part2)
  
  library[[i]] <- readr::read_tsv(lib_file_name,
    col_types = readr::cols(.default = "?", 
                            sgRNA1_Chr = "c", 
                            sgRNA2_Chr = "c", 
                            sgRNA1_WGE_ID = "c", 
                            sgRNA2_WGE_ID = "c"), 
    show_col_types = FALSE)
  
  # filter library per values in data
  library[[i]] <-  library[[i]][match(data[[i]]$sgrna, library[[i]]$ID),]
  
  # rename columns
  data[[i]] <- data[[i]] %>%
    dplyr::rename(ID = sgrna) %>%
    dplyr::mutate(lib = curr_lib, .after = ID) %>%
    dplyr::mutate(ID_lib = paste(ID, curr_lib, sep = "_"), .after = lib) %>%
    dplyr::select(-dplyr::ends_with("_logFC"))
  
  tmp <-  library[[i]] %>%
    dplyr::select(ID, sgRNA1_WGE_ID, sgRNA1_WGE_Sequence, sgRNA2_WGE_ID, sgRNA2_WGE_Sequence) %>%
    dplyr::mutate(SEQ_pair = paste0(sgRNA1_WGE_Sequence, "~", sgRNA2_WGE_Sequence))
  
  data[[i]] <- left_join(data[[i]], tmp, by = "ID") %>%
    dplyr::filter(!duplicated(SEQ_pair))
  
}
names(data) <- names(library) <- lib_name

# which cell lines?
sample_names <- lapply(data, function(x) 
  colnames(x)[grepl("_LFC_RAW", colnames(x))])

sample_names <- lapply(sample_names, function(x) 
  unique(str_split_fixed(string = x, pattern = "[_]", n = 4)[,1]))

df <- data.frame(model_id_CMP = unlist(sample_names), 
                 lib = unname(unlist(mapply(function(x,y) rep(x,length(y)), 
                                            x = lib_name, y = sample_names, SIMPLIFY = T))))

CMP_table <- read_csv("https://cog.sanger.ac.uk/cmp/download/model_list_20230505.csv") %>%
  dplyr::select(model_id, sample_id, model_name, synonyms, tissue, cancer_type, 
                tissue_status, COSMIC_ID, BROAD_ID, CCLE_ID) %>%
  dplyr::rename(model_id_CMP = model_id, 
                sample_id_CMP = sample_id, 
                model_name_CMP = model_name)

model_encore_table <- left_join(df, CMP_table, by = "model_id_CMP") %>%
  mutate(model_name_uppercase = str_replace_all(model_name_CMP, "[-]", "")) %>%
  mutate(model_name_uppercase = toupper(model_name_uppercase))

model_encore_table <- model_encore_table %>% 
  arrange(lib)
CL <- unique(model_encore_table$model_name_CMP)
tab_count <- table(model_encore_table$model_name_CMP, model_encore_table$lib)
tab_count <- tab_count[match(CL, rownames(tab_count)),]
CL_summary <- pheatmap::pheatmap(tab_count,
                                 color = c("grey80", "black"),
                                 breaks = c(0, 0.5, 1),
                                 cluster_rows = TRUE,
                                 cluster_cols = FALSE,
                                 treeheight_row = 0)

### NOT PROPER BATCH CORRECTION! ###
# batch correction and PCA only for common pairs
# plot PCA (opposite size, impossible to correct all)
common_COLO_os <- pca_commonpairs_oppositesize_function(data[1:3]) 
common_BRCA_os <- pca_commonpairs_oppositesize_function(data[4:6])

# get all corrected dataset
data_COLO <- adjust_alldata(list_df = data[1:3])
data_BRCA <- adjust_alldata(list_df = data[4:6])

# plot distribution
plot_CL_distribution(original = data_COLO$original, 
                     adjusted = data_COLO$adj, 
                     common_pairs = data_COLO$combat$common_pairs)

plot_CL_distribution(original = data_BRCA$original, 
                     adjusted = data_BRCA$adj, 
                     common_pairs = data_BRCA$combat$common_pairs)


# get correlation
plot_dist_commonpairs(list_df = data[1:3], oppositesize = TRUE)
plot_dist_commonpairs(list_df = data[4:6], oppositesize = TRUE)

### does it make sense to also correct per CL?
### so far not working...
data_COLO_percl <- adjust_alldata_percl(list_df = data[1:3], 
                                        mat_harm = data_COLO$adj)
plot_CL_distribution(original = data_COLO_percl$original, 
                     adjusted = data_COLO_percl$adj, 
                     common_pairs = data_COLO_percl$combat$common_pairs)

### PROPER BATCH CORRECTION! ###
# batch correction
# plot PCA
common_COLO <- pca_commonpairs_function(data[1:3])
plot_dist_commonpairs(list_df = data[1:3], oppositesize = FALSE)

common_BRCA <- pca_commonpairs_function(data[4:6])
plot_dist_commonpairs(list_df = data[4:6], oppositesize = FALSE)

# correct per PCs
common_COLO_PCcorr <- pca_commonpairs_function(data[1:3], 
                                               correct_for_PCs = TRUE)
plot_dist_commonpairs(mat_raw_corrected = list(raw = t(common_COLO_PCcorr$raw_PCcorr), 
                                               corrected = t(common_COLO_PCcorr$adjusted_PCcorr)))

common_BRCA_PCcorr <- pca_commonpairs_function(data[4:6], 
                                               correct_for_PCs = TRUE)
plot_dist_commonpairs(mat_raw_corrected = list(raw = t(common_BRCA_PCcorr$raw_PCcorr), 
                                               corrected = t(common_BRCA_PCcorr$adjusted_PCcorr)))

# get all corrected dataset
data_COLO <- adjust_alldata(list_df = data[1:3])
data_BRCA <- adjust_alldata(list_df = data[4:6])

# plot distribution
plot_CL_distribution(original = data_COLO$original, 
                     adjusted = data_COLO$adj, 
                     common_pairs = data_COLO$combat$common_pairs)

plot_CL_distribution(original = data_BRCA$original, 
                     adjusted = data_BRCA$adj, 
                     common_pairs = data_BRCA$combat$common_pairs)


