# batch integration
library(tidyverse)
library(sva)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(pROC)
library(CRISPRcleanR)
source("R/batch_integration_functions.R")

#### 
fold <- "/group/iorio/lucia/datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/"
fold_cmp <- "/group/iorio/lucia/datasets/CMP_PORTAL/"
fold_input <- sprintf("%sORIGINAL/", fold)
fold_lib <- sprintf("%sLIBS/", fold)
fold_output <- sprintf("%sBATCH_CORRECTED/", fold)
####

lib_name <- c("COLO1", "COLO2", "COLO3", 
              "BRCA1", "BRCA2", "BRCA3")
data <- library <- list()

for (i in 1:length(lib_name)) {
  
  curr_lib <- lib_name[i]
  print(curr_lib)
  data[[i]] <- read_csv(sprintf("%s%s_FINAL_EXACT_logFC_sgRNA_scaled.txt", fold_input, curr_lib))  
  
  # load lib
  if (grepl("COLO", curr_lib)) {
    part1 <- "COREAD"
  }else{
    part1 <- "BRCA"
  }
  part2 <- stringr::str_sub(curr_lib, start = 5, end = 5)
  lib_file_name <- sprintf("%sENCORE_GI_%s_Library_%s.txt", fold_lib, part1, part2)
  
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
  
  # get same order as in data (possible removal of the same SEQ_pair)
  library[[i]] <-  library[[i]][match(data[[i]]$ID, library[[i]]$ID),] %>%
    dplyr::mutate(lib = curr_lib, .after = ID) %>%
    dplyr::mutate(ID_lib = paste(ID, curr_lib, sep = "_"), .after = lib) %>%
    dplyr::mutate(SEQ_pair = paste0(sgRNA1_WGE_Sequence, "~", sgRNA2_WGE_Sequence))
  
  if ("MyNotes" %in% colnames(library[[i]])) {
    library[[i]] <- library[[i]] %>%
      dplyr::select(-MyNotes)
  }
    
  
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

CMP_table <- read_csv(sprintf("%smodel_annotation/model_list_20230801.csv", fold_cmp)) %>% 
  dplyr::select(model_id, sample_id, model_name, synonyms, tissue, cancer_type, 
                tissue_status, COSMIC_ID, BROAD_ID, CCLE_ID) %>%
  dplyr::rename(model_id_CMP = model_id, 
                sample_id_CMP = sample_id, 
                model_name_CMP = model_name)

model_encore_table <- dplyr::left_join(df, CMP_table, by = "model_id_CMP") %>%
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

# plot PCA:
common_COLO <- pca_commonpairs_function(data[1:3], 
                                        outfold = sprintf("%sCOLO_", fold_output), 
                                        save_plot = TRUE, 
                                        show_plot = TRUE) 

common_BRCA <- pca_commonpairs_function(data[4:6], 
                                        outfold = sprintf("%sBRCA_", fold_output), 
                                        save_plot = TRUE, 
                                        show_plot = TRUE) 

plot_dist_commonpairs(list_df = data[1:3], 
                      outfold = sprintf("%sCOLO_", fold_output))
plot_dist_PPV(list_df = data[1:3], 
              outfold = sprintf("%sCOLO_", fold_output))

plot_dist_commonpairs(list_df = data[4:6], 
                      outfold = sprintf("%sBRCA_", fold_output))
plot_dist_PPV(list_df = data[4:6], 
              outfold = sprintf("%sBRCA_", fold_output))

#### Test Neighbors strategy ####
# save plots
validation_COLO <- validate_NN_approximation(list_df = data[1:3], 
                                             outfold = sprintf("%sCOLO_", fold_output))

validation_BRCA <- validate_NN_approximation(list_df = data[4:6], 
                                             outfold = sprintf("%sBRCA_", fold_output))

# get all corrected dataset
data_COLO <- adjust_alldata_kNN(list_df = data[1:3], 
                                kNN = 5,
                                outfold = sprintf("%sCOLO_", fold_output), 
                                save_plot = TRUE, 
                                show_plot = TRUE) 

data_BRCA <- adjust_alldata_kNN(list_df = data[4:6], 
                                kNN = 5, 
                                outfold = sprintf("%sBRCA_", fold_output), 
                                save_plot = TRUE, 
                                show_plot = TRUE) 

# plot distribution
COLO_allCLs <- plot_CL_distribution(original = data_COLO$original, 
                                    adjusted = data_COLO$adj, 
                                    common_pairs = data_COLO$combat$common_pairs, 
                                    outfold = sprintf("%sCOLO_", fold_output), 
                                    save_plot = TRUE, 
                                    show_plot = TRUE) 

BRCA_allCLs <- plot_CL_distribution(original = data_BRCA$original, 
                                    adjusted = data_BRCA$adj, 
                                    common_pairs = data_BRCA$combat$common_pairs, 
                                    outfold = sprintf("%sBRCA_", fold_output), 
                                    save_plot = TRUE, 
                                    show_plot = TRUE) 

# create final tables and save
data_adj_COLO <- get_complete_table(
  list_df = data[1:3], 
  list_matrix = data_COLO$adj
)
data_or_COLO <- get_complete_table(
  list_df = data[1:3], 
  list_matrix = data_COLO$original
)

# save adjusted output
write.table(file = sprintf("%sCOLO_FINAL_EXACT_logFC_sgRNA_ComBatCorrectionLIBs.txt", fold_output), 
            x = data_adj_COLO, 
            quote = F, 
            col.names = T, 
            row.names = F, 
            sep = "\t")

# save combined libraries 
write.table(file = sprintf("%sENCORE_GI_COREAD_Library_ALL.txt", fold_output), 
            x = bind_rows(library[1:3]), 
            quote = F, 
            col.names = T, 
            row.names = F, 
            sep = "\t")

# save parameters for each guide pair
param <- data_COLO$param_all
save(param, 
     file = sprintf("%sCOLO_FINAL_EXACT_logFC_sgRNA_ComBatParam.RData", fold_output))

data_adj_BRCA <- get_complete_table(
  list_df = data[4:6], 
  list_matrix = data_BRCA$adj
)

data_or_BRCA <- get_complete_table(
  list_df = data[4:6], 
  list_matrix = data_BRCA$original
)

# save adjusted output
write.table(file = sprintf("%sBRCA_FINAL_EXACT_logFC_sgRNA_ComBatCorrectionLIBs.txt", fold_output), 
            x = data_adj_BRCA, 
            quote = F, 
            col.names = T, 
            row.names = F, 
            sep = "\t")

# save combined libraries 
write.table(file = sprintf("%sENCORE_GI_BRCA_Library_ALL.txt", fold_output), 
            x = bind_rows(library[4:6]), 
            quote = F, 
            col.names = T, 
            row.names = F, 
            sep = "\t")

param <- data_BRCA$param_all
save(param, 
     file = sprintf("%sBRCA_FINAL_EXACT_logFC_sgRNA_ComBatParam.RData", fold_output))

### external validation: check CL specific distributions before and after correction ###
ktest_COLO <- test_distributions_per_class(data_adj = data_adj_COLO, 
                                           data_or = data_or_COLO, 
                                           outfold = sprintf("%sCOLO_", fold_output), 
                                           save_plot = TRUE, 
                                           show_plot = TRUE)

ktest_BRCA <- test_distributions_per_class(data_adj = data_adj_BRCA, 
                                           data_or = data_or_BRCA, 
                                           outfold = sprintf("%sBRCA_", fold_output), 
                                           save_plot = TRUE, 
                                           show_plot = TRUE)


# those in library singletons, do they have an in balance in essential genes?
data("ADaM2021_essential")
COLO_lib_genes <- plot_library_genes(data_adj = data_adj_COLO, 
                   data_or = data_or_COLO, 
                   essential_genes = ADaM2021_essential,  
                   outfold = sprintf("%sCOLO_", fold_output), 
                   save_plot = TRUE, 
                   show_plot = TRUE)

BRCA_lib_genes <- plot_library_genes(data_adj = data_adj_BRCA, 
                   data_or = data_or_BRCA, 
                   essential_genes = ADaM2021_essential,  
                   outfold = sprintf("%sBRCA_", fold_output), 
                   save_plot = TRUE, 
                   show_plot = TRUE)







