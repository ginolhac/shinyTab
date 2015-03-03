

# read data from individual data frame
read_data <- function(times, conditions, replicates){
  
  df <- c()
  for(t in times){
    for(cond in conditions){
      for(r in replicates){
        cat(paste("set", t, "h", cond, "#", r, "\n"))
        # use get() to retrieve the object from its name
        # and convert to the Wickham's data frame
        dat <- get(paste("d", t, "_", cond, "_", r, sep = "")) %>% 
          tbl_df() %>%
          select(Accession, Score, Coverage) %>%
          mutate(Score = as.numeric(Score),
                 Coverage = as.numeric(Coverage))
        tmp <- cbind(dat, Condition = rep(cond, nrow(dat)),
                     Time = rep(t, nrow(dat)),
                     Replicate = rep(r, nrow(dat)))
        df <- tbl_df(rbind(df, tmp))
      }
    }
  }
  return(df)
}


plot.PCA <- function(mat){
  
  pca <- prcomp(t(mat), scale = TRUE)
  cond <- colnames(mat)
  pca.df <- data.frame(pca$x[, 1:2], cond = cond)
  pca.df <- cbind(pca.df, colsplit(pca.df$cond, "_", c("condition", "replicate")))
  # get % variances explained
  var_pca <- cumsum((pca$sdev)^2) / sum(pca$sdev^2)
  pp <- ggplot(pca.df)+
    geom_point(aes(x = PC1, y = PC2, colour = condition), alpha = 0.8, size = 5)+
    geom_text(aes(x = PC1 + 0.5, y = PC2 + 0.5, label = replicate))+
    xlab(paste("PC1 (", round(var_pca[1] * 100, 2)," %)"))+
    ylab(paste("PC2 (", round((var_pca[2] - var_pca[1]) * 100, 2)," %)"))+
    ggtitle("principal component analysis")
  return(pp)
}

plot.mds <- function(mat){
  mdsdat <- cmdscale(dist(t(mat)), k = 2)
  mdsdat <- as.data.frame(mdsdat) %>% add_rownames
  mdsdat <- cbind(mdsdat, colsplit(mdsdat$rowname, "_", c("condition", "replicate")))
  p <- ggplot(data = mdsdat)+
    geom_point(aes(x = V1, y = V2, colour = condition), alpha = 0.8, size = 5)+
    geom_text(aes(x = V1 + 0.2, y = V2 + 0.2, label = replicate))+
    ggtitle("multidimensional scaling analysis")
  return(p)
}

# code from Enrico Glaab
contrasts_limma <- function(dat, groups, ID, control, treatment){

  design <- model.matrix(~ -1 + 
                        factor(letters[ifelse(groups[c(which(groups == control),
                               which(groups == treatment))] == control, 0, 1) + 1]))
  
  outnum <- as.vector(ifelse(groups[c(which(groups == control),
                                      which(groups == treatment))] == control,
                             0, 1))
  colnames(design) <- unique(letters[outnum + 1])
  
  corfit <- duplicateCorrelation(dat[,c(which(groups == control),
                                        which(groups == treatment))],
                                 design,
                                 block = ID)
  #corfit$consensus
  fit <- lmFit(dat[,c(which(groups == control),
                         which(groups == treatment))],
               design,
               block = ID,
               correlation = corfit$consensus)
  
  cm <- makeContrasts(comp = b - a, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  return(eBayes(fit2))
}

# wrap the topTAble function to get accession names
topAccession <- function(fit, dat, nb){
  df <- topTable(fit, n = nb) %>% 
    add_rownames()
  accession <- dat %>% 
    add_rownames() %>%
    select(rowname, Accession)
  left_join(df, accession) %>%
    select(Accession, logFC:B) %>%
    return()
}
