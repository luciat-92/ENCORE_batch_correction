# batch integration
library(tidyverse)
library(sva)
library(CRISPRcleanRatSquared)
setwd("/group/iorio/lucia/datasets/ENCORE_SAMPLES_COPYNUMBER/DATA_FREEZE_v4/")

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
  
  data[[i]] <- left_join(data[[i]], tmp, by = "ID")
  
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
# modify combat function
ComBatCP <- function(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                      mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam"),
                     empBayes=TRUE) {
  ## make batch a factor and make a set of indicators for batch
  if(mean.only==TRUE){
    message("Using the 'mean only' version of ComBat")
  }
  if(length(dim(batch))>1){
    stop("This version of ComBat only allows one batch variable")
  }  ## to be updated soon!
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  if (!is.null(ref.batch)){
    ## check for reference batch, check value, and make appropriate changes
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not one of the levels of the batch variable")
    }
    cat("Using batch =",ref.batch, "as a reference batch (this batch won't change)\n")
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  message("Found", nlevels(batch), "batches")
  
  ## A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  } # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){
    mean.only=TRUE
    message("Note: one batch has only one sample, setting mean.only=TRUE")
  }
  n.array <- sum(n.batches)
  ## combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  ## check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  if(!is.null(ref)){
    check[ref] <- FALSE
  } ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  
  ## Number of covariates or covariate levels
  message("Adjusting for", ncol(design)-ncol(batchmod), 'covariate(s) or covariate level(s)')
  
  ## Check if the design is confounded
  if(qr(design)$rank < ncol(design)) {
    ## if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)) {
      stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
    }
    if(ncol(design)>(n.batch+1)) {
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      } else {
        stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
      }
    }
  }
  
  ## Check for missing values
  NAs <- any(is.na(dat))
  if(NAs){
    message(c('Found',sum(is.na(dat)),'Missing Data Values'), sep=' ')}
  ## print(dat[1:2,])
  
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){
    B.hat <- solve(crossprod(design), tcrossprod(t(design), as.matrix(dat)))
  } else { 
    B.hat <- apply(dat, 1, Beta.NA, design) # FIXME
  }
  
  ## change grand.mean for ref batch
  if(!is.null(ref.batch)){
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  
  ## change var.pooled for ref batch
  if (!NAs){
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- ((ref.dat-t(design[batches[[ref]], ] %*% B.hat))^2) %*% rep(1/n.batches[ref],n.batches[ref]) # FIXME
    } else {
      var.pooled <- ((dat-t(design %*% B.hat))^2) %*% rep(1/n.array,n.array) # FIXME
    }
  } else {
    if(!is.null(ref.batch)) {
      ref.dat <- dat[, batches[[ref]]]
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)
    } else {
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)
    }
  }
  
  stand.mean <- t(grand.mean) %*% t(rep(1,n.array)) # FIXME
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp %*% B.hat) #FIXME
  }  
  s.data <- (dat-stand.mean)/(sqrt(var.pooled) %*% t(rep(1,n.array))) # FIXME
  
  ##Get regression batch effect parameters
  message("Fitting L/S model and finding priors")
  batch.design <- design[, 1:n.batch]
  if (!NAs){
    gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                           as.matrix(s.data)))
  } else{
    gamma.hat <- apply(s.data, 1, Beta.NA, batch.design) # FIXME
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE) {
      delta.hat <- rbind(delta.hat,rep(1,nrow(s.data))) 
    } else {
      delta.hat <- rbind(delta.hat, rowVars(s.data[,i], na.rm=TRUE))
    }
  }
  
  if(empBayes){
    ##Find Priors
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apply(delta.hat, 1, sva:::aprior) # FIXME 
    b.prior <- apply(delta.hat, 1, sva:::bprior) # FIXME
    
    ## Plot empirical and parametric priors
    
    if (prior.plots && par.prior) {
      par(mfrow=c(2,2))
      
      ## Top left
      tmp <- density(gamma.hat[1,])
      plot(tmp,  type='l', main=expression(paste("Density Plot of First Batch ",  hat(gamma))))
      xx <- seq(min(tmp$x), max(tmp$x), length=100)
      lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
      
      ## Top Right
      qqnorm(gamma.hat[1,], main=expression(paste("Normal Q-Q Plot of First Batch ", hat(gamma))))
      qqline(gamma.hat[1,], col=2)
      
      ## Bottom Left
      tmp <- density(delta.hat[1,])
      xx <- seq(min(tmp$x), max(tmp$x), length=100)
      tmp1 <- list(x=xx, y=dinvgamma(xx, a.prior[1], b.prior[1]))
      plot(tmp, typ="l", ylim=c(0, max(tmp$y, tmp1$y)),
           main=expression(paste("Density Plot of First Batch ", hat(delta))))
      lines(tmp1, col=2)
      
      ## Bottom Right
      invgam <- 1/qgamma(1-ppoints(ncol(delta.hat)), a.prior[1], b.prior[1])
      qqplot(invgam, delta.hat[1,],
             main=expression(paste("Inverse Gamma Q-Q Plot of First Batch ", hat(delta))),
             ylab="Sample Quantiles", xlab="Theoretical Quantiles")
      lines(c(0, max(invgam)), c(0, max(invgam)), col=2)
    }
    
    ## Find EB batch adjustments
    
    gamma.star <- delta.star <- matrix(NA, nrow=n.batch, ncol=nrow(s.data))
    if (par.prior) {
      message("Finding parametric adjustments")
      results <- bplapply(1:n.batch, function(i) {
        if (mean.only) {
          gamma.star <- postmean(gamma.hat[i,], gamma.bar[i], 1, 1, t2[i])
          delta.star <- rep(1, nrow(s.data))
        }
        else {
          temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i, ],
                               delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                               b.prior[i])
          gamma.star <- temp[1, ]
          delta.star <- temp[2, ]
        }
        list(gamma.star=gamma.star, delta.star=delta.star)
      }, BPPARAM = BPPARAM)
      for (i in 1:n.batch) {
        gamma.star[i,] <- results[[i]]$gamma.star
        delta.star[i,] <- results[[i]]$delta.star
      }
    }
    else {
      message("Finding nonparametric adjustments")
      results <- bplapply(1:n.batch, function(i) {
        if (mean.only) {
          delta.hat[i, ] = 1
        }
        temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                           gamma.hat[i, ], delta.hat[i, ])
        list(gamma.star=temp[1,], delta.star=temp[2,])
      }, BPPARAM = BPPARAM)
      for (i in 1:n.batch) {
        gamma.star[i,] <- results[[i]]$gamma.star
        delta.star[i,] <- results[[i]]$delta.star
      }
    }
  }else{
    #no empirical bayes adjustment:
    gamma.star<-gamma.hat
    delta.star<-delta.hat
  }
  if(!is.null(ref.batch)){
    gamma.star[ref,] <- 0  ## set reference batch mean equal to 0
    delta.star[ref,] <- 1  ## set reference batch variance equal to 1
  }
  
  ## Normalize the Data ###
  message("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j]))) # FIXME
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean # FIXME
  
  ## tiny change still exist when tested on bladder data
  ## total sum of change within each batch around 1e-15 
  ## (could be computational system error).  
  ## Do not change ref batch at all in reference version
  if(!is.null(ref.batch)){
    bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
  }
  
  return(list(correctedData = bayesdata,
              batchDesign = batch.design,
              gamma.star = gamma.star,
              delta.star = delta.star,
              varpool = var.pooled,
              stdmean = stand.mean))
}

combat_correction <- function(list_df){
  
  common_pairs <- base::Reduce(intersect, lapply(list_df, function(x) unique(x$SEQ_pair)))
  data_common <- lapply(list_df, function(x) 
    x[match(common_pairs , x$SEQ_pair), grepl("_LFC_RAW",colnames(x))])
  mat <- list()
  for (i in 1:length(data_common)) {
    
    tmp <- data_common[[i]]
    colnames(tmp) <- str_split_fixed(string = colnames(tmp), pattern = "[_]", n = 4)[,1]
    colnames(tmp) <- paste0(
      model_encore_table$model_name_CMP[match(colnames(tmp), model_encore_table$model_id_CMP)], 
      "_", names(data_common)[i])
    mat[[i]] <- tmp 
  }
  
  tot_mat <- do.call(cbind, mat)
  rownames(tot_mat) <- common_pairs
  ComBat_res <- ComBatCP(dat = as.matrix(tot_mat), 
                      batch = str_split_fixed(colnames(tot_mat), pattern = "_", n = 2)[,2])
  corrected <- as.data.frame(ComBat_res$correctedData)
  
  return(list(raw = tot_mat, 
              corrected = corrected, 
              ComBat_res = ComBat_res))
}

# distribution same pairs
pca_common_function <- function(list_df){
  
  res <- combat_correction(list_df)  
  corrected_common <- res$corrected
  raw_common <- res$raw
  
  # pc from corrected data
  pca_common <- prcomp(t(corrected_common))
  pc_corr <- data.frame(pca_common$x[, 1:4], 
                   CL = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,1],
                   lib = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,2])
  pc_corr$lib <- factor(pc_corr$lib)
  pc_corr$CL <- factor(pc_corr$CL)
  
  pl1 <- ggplot(pc_corr, aes(x = PC1, 
                             y = PC2, 
                             col = CL, 
                             label = CL, 
                             shape = lib)) + 
    geom_point(size = 2) +
    geom_text_repel(size = 2) +
    guides(color = "none") +
    theme_bw() + 
    ggtitle("ComBat corrected", 
            subtitle = paste0("N. common guide pairs: ", length(common_pairs)))
  
  # pc from raw data
  pca_common <- prcomp(t(raw_common))
  pc_raw <- data.frame(pca_common$x[, 1:4], 
                   CL = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,1],
                   lib = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,2])
  pc_raw$lib <- factor(pc_raw$lib)
  pc_raw$CL <- factor(pc_raw$CL)
  pl2 <- ggplot(pc_raw, aes(x = PC1, 
                            y = PC2, 
                            col = CL, 
                            label = CL, 
                            shape = lib)) + 
    geom_point(size = 2) +
    geom_text_repel(size = 2) +
    theme_bw() + 
    guides(color = "none") +
    ggtitle("Raw logFC", 
            subtitle = paste0("N. common guide pairs: ", length(common_pairs)))
  pl <- ggpubr::ggarrange(plotlist = list(pl2, pl1), ncol = 2, common.legend = TRUE)
  print(pl)
  
  return(list(common_pairs = rownames(combat_correction),
              pc_corr = pc_corr,
              pc_raw = pc_raw))
  
}

common_COLO <- pca_common_function(data[1:3])
common_BRCA <- pca_common_function(data[4:6])

BatchCorrection <- function(list_df,
                            stdPrior=TRUE){
  
  
  
  site <- str_split_fixed(colnames(tot_mat), pattern = "_", n = 2)[,2]
  site <- unlist(lapply(1:length(list_df), function(x) rep(nrow(list_df[x])))
    c(rep(site1,ncol(data1)),rep(site2,ncol(data2)))
  adjusted <- AdjustNewData(list_df,combat_res,site,stdPrior)
  return(adjusted)
}

### TODO: emulate BatchCorrection function 
### in https://github.com/DepMap-Analytics/IntegratedCRISPR/blob/main/Combat_HKfunctions.R




# plot distribution for each cell line
plot_CL_distribution <- function(list_df, common_pairs){
  
  tmp_model <- model_encore_table[model_encore_table$lib %in% names(list_df),]
  CL_names <- unique(tmp_model$model_name_CMP)
  pl <- list()
  for (i in 1:length(CL_names)) {
    
    CL_name <- CL_names[i]
    CL_CMP <- unique(tmp_model$model_id_CMP[tmp_model$model_name_CMP == CL_name])
    print(sprintf("%s (%s)", CL_name, CL_CMP))
    
    df_tmp <- lapply(list_df, function(x) 
      x[,c(1:6,grep(CL_CMP, colnames(x)), grep("SEQ_pair", colnames(x)))])
    
    id_keep <- which(sapply(df_tmp, ncol) == 8)
    df_tmp <- df_tmp[id_keep] 
    df_tmp <- lapply(df_tmp,  function(x) 
      x %>% dplyr::rename_at(7, ~"logFC")) 
    df_CL <- bind_rows(df_tmp)
    print(dim(df_CL))
    
    
    pl1 <- ggplot(subset(df_CL, SEQ_pair %in% common_pairs), 
                  aes(x = Note1, y = logFC, fill = lib)) + 
      geom_boxplot(outlier.size = 1) + 
      theme_bw() + 
      xlab("") +
      coord_flip() + 
      ggtitle(sprintf("%s: only common pairs", CL_name))
    
    pl2 <- ggplot(df_CL,
                  aes(x = Note1, y = logFC, fill = lib)) + 
      geom_boxplot(outlier.size = 1) + 
      theme_bw() + 
      xlab("") +
      coord_flip() + 
      ggtitle(sprintf("%s: all pairs", CL_name))
    pl[[i]] <- ggarrange(plotlist = list(pl1, pl2), ncol = 2, common.legend = TRUE)
    print(pl[[i]])
  }
  
}


plot_CL_distribution(list_df = data[1:3], common_pairs = common_COLO$common_pairs)
plot_CL_distribution(list_df = data[4:6], common_pairs = common_BRCA$common_pairs)

