# modify combat function from DepMap -  Sanger integration paper
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
      # plot(tmp, typ="l", ylim = c(0, max(tmp$y, tmp1$y)),
      plot(tmp, typ="l",
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

# get common CLs across a list of screens
# extract only raw logFCs and convert to matrices
harmonize_per_CL <- function(list_df){
  
  # rename columns
  mat <- list()
  for (i in 1:length(list_df)) {
    tmp <- as.data.frame(list_df[[i]])
    rownames(tmp) <- tmp$SEQ_pair
    tmp <- tmp[, grepl("_LFC_RAW",colnames(tmp))]
    colnames(tmp) <- str_split_fixed(string = colnames(tmp), pattern = "[_]", n = 4)[,1]
    colnames(tmp) <- model_encore_table$model_name_CMP[match(colnames(tmp), model_encore_table$model_id_CMP)]
    mat[[i]] <- t(tmp) 
  }
  
  common_CL <- base::Reduce(intersect,lapply(mat, rownames))
  data_common <- lapply(mat, function(x) x[common_CL, ])
  
  names(data_common) <- names(list_df)
  return(data_common)
  
}

# perform combat correction (per guide pairs)
combat_correction <- function(list_df){
  
  list_df_h <- harmonize_per_CL(list_df)
  common_pairs <- base::Reduce(intersect, 
                               lapply(list_df, function(x) unique(x$SEQ_pair)))
  
  data_common <- lapply(list_df_h, function(x) t(x[, common_pairs]))
  for (i in 1:length(data_common)) {
    colnames(data_common[[i]]) <- paste(colnames(data_common[[i]]), names(list_df)[i], sep = "_")
  }
  
  tot_mat <- do.call(cbind, data_common)
  ComBat_res <- ComBatCP(
    dat = as.matrix(tot_mat),
    batch = str_split_fixed(colnames(tot_mat), pattern = "_", n = 2)[,2])
  corrected <- as.data.frame(ComBat_res$correctedData)

  # annotate combat res and plot
  df_annot <- mapply(function(x, y) 
    x[match(rownames(y), x$SEQ_pair), !grepl("LFC_RAW", colnames(x))], 
    x = data[names(data_common)], y = data_common, SIMPLIFY = FALSE)
  df_annot_unique <- data.frame(SEQ_pair = df_annot[[1]]$SEQ_pair, 
                                Gene_Pair = df_annot[[1]]$Gene, 
                                Note1 = apply(sapply(df_annot, function(x) x$Note1), 
                                              1, function(y) paste0(sort(unique(y)), collapse = ",")), 
                                Note2 = apply(sapply(df_annot, function(x) x$Note2), 
                                              1, function(y) paste0(sort(unique(y)), collapse = ",")))
  
  colnames(ComBat_res$gamma.star) <- colnames(ComBat_res$delta.star) <- df_annot_unique$SEQ_pair
  rownames(ComBat_res$gamma.star) <- rownames(ComBat_res$delta.star)  <- names(df_annot)
  
  df_ComBat_param_gamma <- melt(ComBat_res$gamma.star, 
                          varnames = c("lib", "SEQ_pair")) %>%
    dplyr::mutate(param = "gamma")
  
  df_ComBat_param_delta <- melt(ComBat_res$delta.star, 
                                varnames = c("lib", "SEQ_pair")) %>%
    dplyr::mutate(param = "delta") 
  df_ComBat_param <- rbind(df_ComBat_param_gamma, df_ComBat_param_delta) %>%
    left_join(df_annot_unique, by = "SEQ_pair")
  
  # plot dist of parameters
  pl_lib <- ggplot(df_ComBat_param, 
                aes(x = lib, y = value, fill = lib)) + 
    geom_violin() + 
    geom_boxplot(fill = "white", outlier.size = 1, width = 0.2) + 
    theme_bw() + 
    facet_wrap(.~param, scales = "free_y") +
    xlab("") +
    ylab("ComBat param estimates")
  
  pl_class <- ggplot(df_ComBat_param, 
               aes(x = Note1, y = value, fill = lib)) + 
    geom_boxplot(outlier.size = 1) + 
    theme_bw() + 
    facet_wrap(.~param, scales = "free_x") +
    xlab("") + 
    ylab("ComBat param estimates") +
    coord_flip()
  
  print(pl_lib)
  print(pl_class)
  
  return(list(raw = tot_mat,
              corrected = corrected,
              ComBat_res = ComBat_res))
}

# PC plot for common pairs
# plot before and after combat correction
pca_commonpairs_function <- function(list_df){
  
  res <- combat_correction(list_df)
  corrected_common <- res$corrected
  raw_common <- res$raw
  common_pairs <- rownames(raw_common)
  
  # pc from corrected data
  pca_common <- prcomp(t(corrected_common), scale. = TRUE)
  pc_corr <- data.frame(pca_common$x,
                        CL = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,1],
                        lib = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,2])
  pc_corr$lib <- factor(pc_corr$lib)
  pc_corr$CL <- factor(pc_corr$CL)
  
  pl1 <- ggplot(pc_corr, aes(x = PC1,
                             y = PC2,
                             col = CL,
                             label = CL,
                             group = CL)) +
    geom_point(size = 2, aes(shape = lib)) +
    geom_line(alpha = 0.5) + 
    geom_text_repel(size = 2) +
    guides(color = "none") +
    theme_bw() +
    ggtitle("ComBat corrected",
            subtitle = paste0("N. common guide pairs: ", length(common_pairs)))
  
  # pc from raw data
  pca_common <- prcomp(t(raw_common), scale. = TRUE)
  pc_raw <- data.frame(pca_common$x,
                       CL = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,1],
                       lib = str_split_fixed(rownames(pca_common$x), pattern = "_", n = 2)[,2])
  pc_raw$lib <- factor(pc_raw$lib)
  pc_raw$CL <- factor(pc_raw$CL)
  pl2 <- ggplot(pc_raw, aes(x = PC1,
                            y = PC2,
                            col = CL,
                            label = CL,
                            group = CL)) +
    geom_point(size = 2, aes(shape = lib)) +
    geom_line(alpha = 0.5) + 
    geom_text_repel(size = 2) +
    theme_bw() +
    guides(color = "none") +
    ggtitle("Raw logFC",
            subtitle = paste0("N. common guide pairs: ", length(common_pairs)))
  pl <- ggpubr::ggarrange(plotlist = list(pl2, pl1), ncol = 2, common.legend = TRUE)
  print(pl)
  
  return(list(common_pairs = common_pairs,
              pc_corr = pc_corr,
              pc_raw = pc_raw))
  
}


# distance between each element in a matrix
dist2 <- function(X,C) {
  # from SNFtools
  
  ndata <- nrow(X)
  ncentres <- nrow(C)
  
  sumsqX <- rowSums(X^2)
  sumsqC <- rowSums(C^2)
  
  XC <- 2 * (X %*% t(C))
  res <- matrix(rep(sumsqX, times = ncentres), ndata, ncentres) + 
    t(matrix(rep(sumsqC, times = ndata), ncentres, ndata)) - XC
  res[res < 0] <- 0
  
  return(res)
}

# approximate combat param via kNN
param_approx_kNN <- function(
  matrix_data, 
  gamma_param, 
  delta_param){
  
  common_pairs <- names(gamma_param)
  unique_pairs <- setdiff(colnames(matrix_data), common_pairs)
  
  # distance between common pairs and unique pairs only
  dist_logFCs <- sqrt(dist2(t(matrix_data[,common_pairs]),
                            t(matrix_data[,unique_pairs])))
  
  closest_guides <- apply(dist_logFCs, 2, function(x) 
    order(x)[1:kNN], 
    simplify = TRUE)
  
  gamma_approx_unique <- apply(closest_guides, 2, function(x)
    mean(gamma_param[common_pairs[x]])
  )
  
  delta_approx_unique <- apply(closest_guides, 2, function(x)
    mean(delta_param[common_pairs[x]])
  )
  
  return(list(
    gamma = gamma_approx_unique, 
    delta = delta_approx_unique
  ))
  
}

# adjust data, combat estimates per guide pairs
adjust_alldata_kNN <- function(list_df, kNN){

combat_res_all <- combat_correction(list_df = list_df)
combat_res_all$common_pairs <- rownames(combat_res_all$ComBat_res$correctedData)
  
combat_res <- combat_res_all$ComBat_res
list_logFC <- harmonize_per_CL(list_df)
n_guide_common <- nrow(combat_res$correctedData)

# get combat param
batch.design <- combat_res$batchDesign
gamma.star <- combat_res$gamma.star # additive batch effect
delta.star <- combat_res$delta.star # multiplicative batch effect
stand.mean <- combat_res$stdmean
var.pooled <- combat_res$varpool

# standardize each guide pair, Mean and Sd computed across CLs
stand.mean_all <- lapply(list_logFC, function(x) rowMeans(t(x)) %*% t(rep(1, ncol(t(x)))))
var.pooled_all <- lapply(list_logFC, function(x) rowSds(t(x)) %*% t(rep(1, ncol(t(x)))))
# NOTE:
# we cannot use stand.mean and var.pooled (stdPrior) because we need to estimate guide pairs not in common
# but correlation is high for common pairs with this and ComBatCP strategy 
# >0.95 for stand.mean and >0.84 for var.pooled
s.data_all <- mapply(function(x,y,z) (t(x) - y)/z, 
                     x = list_logFC, 
                     y = stand.mean_all, 
                     z = var.pooled_all, 
                     SIMPLIFY = FALSE)

# correct based on bayes estimates
mean_common <- mapply(function(x, y) 
  matrix(gamma.star[x,], 
         nrow = n_guide_common, 
         ncol = nrow(y)), 
  x = 1:nrow(gamma.star), y = list_logFC, SIMPLIFY = FALSE)

var_common <- mapply(function(x, y) 
  sqrt(delta.star[x,]) %*% t(rep(1,nrow(y))), 
  x = 1:nrow(delta.star), y = list_logFC, 
  SIMPLIFY = FALSE)

# for those not in common, use kNN
mean_unique <- var_unique <- list()
mean_all <- var_all <- list()
for (id in 1:length(list_logFC)) {
  
  param_approx <- param_approx_kNN(matrix_data = list_logFC[[id]], 
                   gamma_param = gamma.star[id,],
                   delta_param = delta.star[id,])
  
  mean_unique[[id]] <- matrix(param_approx$gamma, 
                              nrow = length(param_approx$gamma), 
                              ncol = nrow(list_logFC[[id]]))
  
  var_unique[[id]] <- matrix(param_approx$delta, 
                              nrow = length(param_approx$delta), 
                              ncol = nrow(list_logFC[[id]]))
  
  rownames(mean_common[[id]]) <- colnames(gamma.star)
  rownames(var_common[[id]]) <- colnames(delta.star)
  rownames(mean_unique[[id]]) <- names(param_approx$gamma)
  rownames(var_unique[[id]]) <- names(param_approx$delta)
  
  # put back together
  var_all[[id]] <- rbind(var_common[[id]], var_unique[[id]])
  mean_all[[id]] <- rbind(mean_common[[id]], mean_unique[[id]])
  # order properly
  var_all[[id]] <- var_all[[id]][rownames(s.data_all[[id]]),]
  mean_all[[id]] <- mean_all[[id]][rownames(s.data_all[[id]]),]
  
}

adjusted.data_all <- mapply(function(x,y,z) (x - y)/z, 
                            x = s.data_all, 
                            y =  mean_all, 
                            z = var_all, 
                            SIMPLIFY = FALSE)

# add original mean and variance back
adjusted.data_all <- mapply(function(x,y,z) t( (x*z) + y ), 
                            x = adjusted.data_all, 
                            y = stand.mean_all, 
                            z = var.pooled_all, 
                            SIMPLIFY = FALSE)

return(list(adj = adjusted.data_all, 
            original = list_logFC, 
            combat = combat_res_all))

}


# plot distributions per CL
plot_CL_distribution <- function(original, adjusted, common_pairs){
  
  CL_names <- rownames(original[[1]])
  # get annotation
  df_annot <- mapply(function(x, y) 
    x[match(colnames(y), x$SEQ_pair), !grepl("LFC_RAW", colnames(x))], 
    x = data[names(original)], y = original, SIMPLIFY = FALSE)
  
  df_annot_common <- lapply(df_annot, function(x) x[match(common_pairs, x$SEQ_pair),])
  # get unique
  df_annot_common_unique <- data.frame(SEQ_pair = df_annot_common[[1]]$SEQ_pair, 
                                Gene_Pair = df_annot_common[[1]]$Gene, 
                                Note1 = apply(sapply(df_annot_common, function(x) x$Note1), 
                                              1, function(y) paste0(sort(unique(y)), collapse = ",")), 
                                Note2 = apply(sapply(df_annot_common, function(x) x$Note2), 
                                              1, function(y) paste0(sort(unique(y)), collapse = ",")))
  
  ## plot everything  
  pl <- list()
  df_CL <- list()
  for (i in 1:length(CL_names)) {
    
    CL_name <- CL_names[i]
    print(sprintf("%s", CL_name))
    
    df_tmp <- mapply(function(x,y,z) 
      cbind(z, data.frame(logFC_or = x[CL_name,, drop = T], logFC_adj = y[CL_name,, drop = T])), 
      x = original, y = adjusted, z = df_annot, SIMPLIFY = F)
    
    df_CL[[i]] <- bind_rows(df_tmp)
    print(dim(df_CL[[i]]))
    
    df_CL[[i]] <- left_join(suffix = c("lib_spec", "lib_common"), 
                            df_CL[[i]], 
                            df_annot_common_unique, 
                            by = "SEQ_pair")
    df_CL[[i]]$Note1 <- df_CL[[i]]$Note1lib_common
    df_CL[[i]]$Note1[is.na(df_CL[[i]]$Note1lib_common)] <- df_CL[[i]]$Note1lib_spec[is.na(df_CL[[i]]$Note1lib_common)]
    
    pl1 <- ggplot(df_CL[[i]], 
                  aes(x = Note1, y = logFC_adj, fill = lib)) + 
      geom_boxplot(outlier.size = 0.5) + 
      theme_bw() + 
      xlab("") +
      coord_flip() + 
      ggtitle(sprintf("%s: adjusted", CL_name))
    
    pl2 <- ggplot(df_CL[[i]],
                  aes(x = Note1, y = logFC_or, fill = lib)) + 
      geom_boxplot(outlier.size = 0.5) + 
      theme_bw() + 
      xlab("") +
      coord_flip() + 
      ggtitle(sprintf("%s: original", CL_name))
    
    pl3 <- ggplot(subset(df_CL[[i]], SEQ_pair %in% common_pairs),
                  aes(x = Note1, y = logFC_adj, fill = lib)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_bw() +
      xlab("") +
      coord_flip() +
      ggtitle(sprintf("%s: adjusted (common pairs)", CL_name))

    pl4 <- ggplot(subset(df_CL[[i]], SEQ_pair %in% common_pairs),
                  aes(x = Note1, y = logFC_or, fill = lib)) +
      geom_boxplot(outlier.size = 0.5) +
      theme_bw() +
      xlab("") +
      coord_flip() +
      ggtitle(sprintf("%s: original (common pairs)", CL_name))

    pl[[i]] <- ggpubr::ggarrange(plotlist = list(pl4, pl3, pl2, pl1), ncol = 2, nrow = 2, align = "hv", common.legend = TRUE)
    print(pl[[i]])
  }
  
  # plot dist across all CLs
  df_tot <- do.call(rbind, df_CL)  
  pl1 <- ggplot(df_tot, 
                aes(x = Note1, y = logFC_adj, fill = lib)) + 
    geom_boxplot(outlier.size = 1) + 
    theme_bw() + 
    xlab("") +
    coord_flip() + 
    ggtitle(sprintf("%s: adjusted", "All CLs"))
  
  pl2 <- ggplot(df_tot,
                aes(x = Note1, y = logFC_or, fill = lib)) + 
    geom_boxplot(outlier.size = 1) + 
    theme_bw() + 
    xlab("") +
    coord_flip() + 
    ggtitle(sprintf("%s: original", "All CLs"))
  
  pl3 <- ggplot(subset(df_tot, SEQ_pair %in% common_pairs), 
                aes(x = Note1lib_common, y = logFC_adj, fill = lib)) + 
    geom_boxplot(outlier.size = 1) + 
    theme_bw() + 
    xlab("") +
    coord_flip() + 
    ggtitle(sprintf("%s: adjusted (common pairs)", "All CLs"))
  
  pl4 <- ggplot(subset(df_tot, SEQ_pair %in% common_pairs),
                aes(x = Note1lib_common, y = logFC_or, fill = lib)) + 
    geom_boxplot(outlier.size = 1) + 
    theme_bw() + 
    xlab("") +
    coord_flip() + 
    ggtitle(sprintf("%s: original (common pairs)", "All CLs"))
  
  pl <- ggpubr::ggarrange(plotlist = list(pl4, pl3, pl2, pl1), ncol = 2, nrow = 2, align = "hv", common.legend = TRUE)
  print(pl)
  
}

# get correlation
dist_commonpairs <- function(mat_common){
  
  dist_mat <- as.matrix(dist(t(mat_common),method = 'euclidean'))
  # cor_res <- cor(mat_common, method = "pearson")
  libs <- str_split_fixed(colnames(dist_mat), pattern = "_", n = 2)[,2]
  libs_name <- unique(libs)
  CLs <- str_split_fixed(colnames(dist_mat), pattern = "_", n = 2)[,1]
  CLs_names <- unique(CLs)
  
  inner_libs <- sapply(libs_name, function(x) median(dist_mat[libs == x, libs == x][upper.tri(dist_mat[libs == x, libs == x])]))
  outer_libs <- sapply(libs_name, function(x) median(dist_mat[libs == x, libs != x]))
  
  inner_CLs <- sapply(CLs_names, function(x) median(dist_mat[CLs == x, CLs == x][upper.tri(dist_mat[CLs == x, CLs == x])]))
  outer_CLs <- sapply(CLs_names, function(x) median(dist_mat[CLs == x, CLs != x]))
  
  return(list(libs = data.frame(name = names(inner_libs), inner = inner_libs, outer = outer_libs), 
              CLs = data.frame(name = names(inner_CLs), inner = inner_CLs, outer = outer_CLs)))
}

plot_dist_commonpairs <- function(list_df){
  
  res_combat <- combat_correction(list_df)
  
  dist_raw <- dist_commonpairs(res_combat$raw)
  dist_combat <- dist_commonpairs(res_combat$corrected)
  
  CL_dist <- data.frame(name = rep(dist_raw$CLs$name, 2), dist = c(dist_raw$CLs$inner, dist_raw$CLs$outer), 
                        type_dist = c(rep("inner", nrow(dist_raw$CLs)), rep("outer", nrow(dist_raw$CLs))), 
                        type_combat = "Raw")
  CL_dist <- rbind(CL_dist, data.frame(name = rep(dist_combat$CLs$name, 2), dist = c(dist_combat$CLs$inner, dist_combat$CLs$outer), 
                                       type_dist = c(rep("inner", nrow(dist_combat$CLs)), rep("outer", nrow(dist_combat$CLs))), 
                                       type_combat = "ComBat corrected"))
  pl1 <- ggplot(CL_dist , aes(x = type_combat,
                              y = dist,
                              fill = type_dist)) +
    geom_boxplot(alpha = 0.4) +
    geom_jitter(position = position_dodge(width = 0.7)) +
    xlab("") + 
    theme_bw() + 
    ggtitle("Distance among and outside same CLs")
  
  libs_dist <- data.frame(name = rep(dist_raw$libs$name, 2), dist = c(dist_raw$libs$inner, dist_raw$libs$outer), 
                          type_dist = c(rep("inner", nrow(dist_raw$libs)), rep("outer", nrow(dist_raw$libs))), 
                          type_combat = "Raw")
  libs_dist <- rbind(libs_dist, data.frame(name = rep(dist_combat$libs$name, 2), dist = c(dist_combat$libs$inner, dist_combat$libs$outer), 
                                           type_dist = c(rep("inner", nrow(dist_combat$libs)), rep("outer", nrow(dist_combat$libs))), 
                                           type_combat = "ComBat corrected"))
  pl2 <- ggplot(libs_dist , aes(x = type_combat,
                                y = dist,
                                fill = type_dist)) +
    geom_boxplot(alpha = 0.4) +
    geom_jitter(position = position_dodge(width = 0.7)) +
    xlab("") + 
    theme_bw() + 
    ggtitle("Distance among and outside same library")
  pl <- ggpubr::ggarrange(plotlist = list(pl2, pl1), ncol = 2, common.legend = TRUE)
  print(pl)
  
  return(list(CL = CL_dist, lib = libs_dist))
  
}

########
# can we estimate gamma and delta from the closest elements?
# validation based on cohen's d inside and outside kNN
get_stat_closest <- function(id_lib, 
                             matrix_data, 
                             combat_param, 
                             kNN = NULL, 
                             quant_prob = NULL){
  
  # explain why euclidean distance makes sense
  dist_guides <- as.matrix(dist(t(matrix_data[[id_lib]])))
  if (is.null(kNN) & is.null(quant_prob)) {
    stop("one between kNN and quant_prob should not be NULL")
  }else{
    # for each guide, find closest guides
    if (!is.null(kNN)) {
      closest_guides <- apply(dist_guides, 1, function(x) order(x)[1:(kNN + 1)], simplify = FALSE)
      closest_guides_dist <- apply(dist_guides, 1, function(x) sort(x)[1:(kNN + 1)], simplify = FALSE)
    }else{
      max_dist <- quantile(as.vector(dist_guides[upper.tri(dist_guides, diag = FALSE)]), 
                           probs = quant_prob)
      closest_guides <- apply(dist_guides, 1, function(x) which(x < max_dist), simplify = FALSE)
      closest_guides_dist <- apply(dist_guides, 1, function(x) x[x < max_dist], simplify = FALSE)
    }
  }
  
  gamma_closest <- lapply(closest_guides, 
                          function(x) combat_param$gamma.star[id_lib,x])
  gamma_compl <- lapply(closest_guides, 
                        function(x) combat_param$gamma.star[id_lib,-x])
  delta_closest <- lapply(closest_guides, 
                          function(x) combat_param$delta.star[id_lib,x])
  delta_compl <- lapply(closest_guides, 
                        function(x) combat_param$delta.star[id_lib,-x])
  
  ###
  n_pairs <- ncol(matrix_data[[id_lib]])
  cohensd_gamma <- sapply(1:n_pairs, function(x) 
    tryCatch(cohens_d(
      x = abs(combat_param$gamma.star[id_lib,x] - gamma_closest[[x]]), 
      y = abs(combat_param$gamma.star[id_lib,x] - gamma_compl[[x]]), 
      pooled_sd = FALSE)$Cohens_d, 
      error = function(e){NA}))
  
  ttest_gamma <- sapply(1:n_pairs, function(x) 
    tryCatch(t.test(
      x = abs(combat_param$gamma.star[id_lib,x] - gamma_closest[[x]]),
      y = abs(combat_param$gamma.star[id_lib,x] - gamma_compl[[x]]), 
      alternative = "less")$p.value, 
      error = function(e){NA}))
  
  cohensd_delta <- sapply(1:n_pairs, function(x) 
    tryCatch(cohens_d(
      x = abs(combat_param$delta.star[id_lib,x] - delta_closest[[x]]), 
      y = abs(combat_param$delta.star[id_lib,x] - delta_compl[[x]]), 
      pooled_sd = FALSE)$Cohens_d, 
      error = function(e){NA}))
  
  ttest_delta <- sapply(1:n_pairs, function(x) 
    tryCatch(t.test(
      x = abs(combat_param$delta.star[id_lib,x] - delta_closest[[x]]), 
      y = abs(combat_param$delta.star[id_lib,x] - delta_compl[[x]]), 
      alternative = "less")$p.value,
      error = function(e){NA}))
  
  df_out <- data.frame(SEQ_pair = rep(colnames(matrix_data[[id_lib]]), 2), 
                       cohens_d = c(cohensd_gamma, cohensd_delta), 
                       ttest_pval = c(ttest_gamma, ttest_delta), 
                       type = c(rep("gamma", n_pairs), rep("delta", n_pairs)), 
                       id_lib = id_lib)
  
  if (!is.null(kNN)) {
    df_out$kNN <- kNN
  }
  
  if (!is.null(quant_prob)) {
    df_out$quant_prob <- quant_prob
    df_out$max_dist <- max_dist
  }
  
  return(df_out)
  
}

# can we estimate gamma and delta from the closest elements?
# validation based on auc using ED distance between param as predictor
compute_auc <- function(id_lib, 
                        matrix_data, 
                        combat_param, 
                        kNN){
  
  # distance of logFCs across CLs
  dist_guides <- as.matrix(dist(t(matrix_data[[id_lib]])))
  closest_guides <- apply(dist_guides, 1, function(x) order(x)[1:(kNN + 1)], simplify = FALSE)
  n_pairs <- ncol(matrix_data[[id_lib]])
  
  # distance of the gamma,delta parameters across guide pairs 
  param_vector <- data.frame(gamma = combat_param$gamma.star[id_lib,], 
                             delta = combat_param$delta.star[id_lib,])
  dist_param <- as.matrix(dist(param_vector, method = "euclidean"))
  
  # response = kNN of the guide pair
  # predictor = distance between the guide pair and all the other guide pairs
  auc_res <- sapply(1:n_pairs, function(x)
    as.numeric(auc(roc(
      response = as.numeric(1:n_pairs %in% closest_guides[[x]]), 
      predictor = unname(dist_param[,x]), 
      quiet = TRUE, 
      direction = ">"))))
  
  df_out <- data.frame(SEQ_pair = colnames(matrix_data[[id_lib]]),
                       auc = auc_res, 
                       id_lib = id_lib, 
                       kNN = kNN)
  return(df_out)
  
}

# combine all validation strategies:
validate_NN_approximation <- function(list_df){
  
  res <- combat_correction(list_df)
  corrected_common <- res$corrected
  raw_common <- res$raw
  common_pairs <- rownames(raw_common)
  matrix_data <- lapply(harmonize_per_CL(list_df), function(x) 
    x[,common_pairs])
  
  # varying kNN
  kNN_list <- c(5, seq(10, 100, by = 10))
  df_kNN <- data.frame()
  for (kNN_p in kNN_list) {
    print(kNN_p)
    tmp <- lapply(1:length(list_df), function(x)
      get_stat_closest(id_lib = x, 
                       matrix_data = matrix_data, 
                       combat_param = res$ComBat_res, 
                       kNN = kNN_p))  
    df_kNN <- rbind(do.call(rbind, tmp), df_kNN)
  }
  
  pl <- ggplot(df_kNN, aes(x = as.factor(kNN),
                     y = cohens_d)) +
    geom_boxplot() +
    facet_wrap(.~type) +
    xlab("k Nearest Neigh.") +
    ylab("Cohen's d:\n |param_i - neigh.| VS |param_i - not neigh.|") +
    theme_bw()
  print(pl)
  
  # varying prob quantile
  quant_list <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  df_quant <- data.frame()
  for (quant_p in quant_list) {
    print(quant_p)
    tmp <- lapply(1:length(list_df), function(x)
      get_stat_closest(id_lib = x, 
                       matrix_data = matrix_data, 
                       combat_param = res$ComBat_res, 
                       quant_prob = quant_p))  
    df_quant  <- rbind(do.call(rbind, tmp), df_quant)
  }
  
  pl <- ggplot(df_quant, aes(x = as.factor(quant_prob),
                       y = cohens_d)) +
    geom_boxplot() +
    facet_wrap(.~type) +
    xlab("Quantile perc. of distance distribution to define maximum Eucl. Distance") +
    ylab("Cohen's d:\n |param_i - neigh.| VS |param_i - not neigh.|") +
    theme_bw()
  print(pl)
  
  count_na <- as.data.frame(table(is.na(df_quant$cohens_d), df_quant$quant_prob)) %>%
    mutate(Var1 = case_when(Var1 == "TRUE" ~ "no guide pairs in the range", 
                            Var1 == "FALSE" ~ "> 0 guide pairs in the range"))
  
  pl <- ggplot(count_na, aes(x = as.factor(Var2),
                       y = Freq, 
                       fill = Var1)) +
    geom_bar(stat = "identity") +
    xlab("Quantile perc. of distance distribution to define maximum Eucl. Distance") +
    ylab("N. guide pairs") + 
    theme_bw() +
    theme(legend.position = "bottom")
  print(pl)
  
  ### Validation 2:
  # create TPR, FPR. Compare the selection of closest kNN with the "closest" param
  auc_kNN <- data.frame()
  for (kNN_p in kNN_list) {
    print(kNN_p)
    tmp <- lapply(1:length(list_df), function(x)
      compute_auc(id_lib = x, 
                  matrix_data = matrix_data, 
                  combat_param = res$ComBat_res, 
                  kNN = kNN_p))  
    auc_kNN <- rbind(do.call(rbind, tmp), auc_kNN)
  }
  
  pl <- ggplot(auc_kNN, aes(x = as.factor(kNN),
                      y = auc)) +
    geom_boxplot() +
    xlab("k Nearest Neigh.") +
    ylab("AUROC") +
    theme_bw()
  print(pl)
  
  return(list(cohensd_kNN = df_kNN, 
              cohensd_range = df_quant, 
              auc_kNN = auc_kNN))
  
}




