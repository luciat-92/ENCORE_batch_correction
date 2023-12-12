# batch integration
library(tidyverse)
library(sva)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(pROC)
library(CRISPRcleanR)
library(effectsize)
source("R/batch_integration_functions.R")

#### 
fold <- "/group/iorio/lucia/datasets/ENCORE_SAMPLES_COPYNUMBER/"
fold_cmp <- "/group/iorio/lucia/datasets/CMP_PORTAL/"
fold_input <- sprintf("%sDATA_FREEZE_v4_NOV_2023/cas9neg_rdsFiles/", fold)
fold_lib <- sprintf("%sDATA_META/LIBS/", fold)
fold_output <- sprintf("%sDATA_FREEZE_v4_NOV_2023/BATCH_CORRECTED/", fold)
####

lib_name <- c("COLO1", "COLO2", "COLO3", 
              "BRCA1", "BRCA2", "BRCA3") %>% tolower()

res_c91 <- load_data_rds(
  fold_input = fold_input, 
  fold_lib = fold_lib, 
  negcontrol = "c91")

res_c92 <- load_data_rds(
  fold_input = fold_input, 
  fold_lib = fold_lib, 
  negcontrol = "c92")

# create avg of the two plasmids
data_cavg <- list()
data_c91 <- res_c91$data
data_c92 <- res_c92$data

s_91 <- lapply(data_c91, function(x) 
  colnames(x)[grepl("SID", colnames(x))])
data_c91_logFC <- mapply(function(x, y) x[, y], x = data_c91, y = s_91)

s_92 <- lapply(data_c92, function(x) 
  colnames(x)[grepl("SID", colnames(x))])
data_c92_logFC <- mapply(function(x, y) x[, y], x = data_c92, y = s_92)

for (i in 1:length(lib_name)) {
  
  lib_idx <- lib_name[i]
  v1 <- data_c91_logFC[[lib_idx]]
  v2 <- data_c92_logFC[[lib_idx]]
  
  if (nrow(v1) != nrow(v2)) {
    stop("Not the same number of guide pairs!")
  }
  if (!identical(data_c91[[i]]$ID, data_c92[[i]]$ID)) {
    stop("Guide pairs do not have the same order!")
  }
  if (!identical(colnames(v1), colnames(v2))) {
    stop("Not same samples in cas9negatives")
  }
  
  data_cavg_logFC <- (v1 + v2)/2
  # use c91, same in both!
  data_cavg[[i]] <- cbind(data_c91[[i]][, !colnames(data_c91[[i]]) %in% s_91[[i]]], data_cavg_logFC)
}
names(data_cavg) <- lib_name
res_cavg <- list(data = data_cavg, library = res_c91$library)

#####################
# which cell lines?
# same for all datasets, use data_cavg
sample_names <- lapply(data_cavg, function(x) 
  colnames(x)[grepl("SID", colnames(x))])

df <- data.frame(model_id_CMP = unlist(sample_names), 
                 lib = unname(unlist(mapply(function(x,y) rep(x,length(y)), 
                                            x = lib_name, 
                                            y = sample_names, SIMPLIFY = T))))

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
############
# plot PCA:
common_COLO <- pca_commonpairs_function(res_cavg$data[1:3], 
                                        # outfold = sprintf("%sCOLO_", fold_output), 
                                        save_plot = FALSE, 
                                        show_plot = TRUE) 

common_BRCA <- pca_commonpairs_function(res_cavg$data[4:6], 
                                        # outfold = sprintf("%sBRCA_", fold_output), 
                                        save_plot = FALSE, 
                                        show_plot = TRUE) 

plot_dist_commonpairs(list_df = res_cavg$data[1:3])
                      # outfold = sprintf("%sCOLO_", fold_output))
plot_dist_PPV(list_df = res_cavg$data[1:3])
              # outfold = sprintf("%sCOLO_", fold_output))

plot_dist_commonpairs(list_df = res_cavg$data[4:6])
                      #outfold = sprintf("%sBRCA_", fold_output))
plot_dist_PPV(list_df = res_cavg$data[4:6])
              # outfold = sprintf("%sBRCA_", fold_output))

#### Test Neighbors strategy ####
# save plots
validation_COLO <- validate_NN_approximation(list_df = res_cavg$data[1:3]) 
                                             # outfold = sprintf("%sCOLO_", fold_output))

validation_BRCA <- validate_NN_approximation(list_df = res_cavg$data[4:6])
                                             # outfold = sprintf("%sBRCA_", fold_output))

# get all corrected dataset
data_COLO <- adjust_alldata_kNN(list_df = res_cavg$data[1:3], 
                                kNN = 5,
                                #outfold = sprintf("%sCOLO_", fold_output), 
                                save_plot = FALSE, 
                                show_plot = TRUE) 

data_BRCA <- adjust_alldata_kNN(list_df = res_cavg$data[4:6], 
                                kNN = 5, 
                                # outfold = sprintf("%sBRCA_", fold_output), 
                                save_plot = FALSE, 
                                show_plot = TRUE) 

# plot distribution
COLO_allCLs <- plot_CL_distribution(original = data_COLO$original, 
                                    adjusted = data_COLO$adj, 
                                    common_pairs = data_COLO$combat$common_pairs, 
                                    # outfold = sprintf("%sCOLO_", fold_output), 
                                    save_plot = FALSE, 
                                    show_plot = TRUE) 

BRCA_allCLs <- plot_CL_distribution(original = data_BRCA$original, 
                                    adjusted = data_BRCA$adj, 
                                    common_pairs = data_BRCA$combat$common_pairs, 
                                    # outfold = sprintf("%sBRCA_", fold_output), 
                                    save_plot = FALSE, 
                                    show_plot = TRUE) 
#### FROM HERE ####
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







