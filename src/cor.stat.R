
filter_by_abun_freq <- function(dt, low_freq = 0.1, low_abun = 1e-04)
{
  freq <- apply(dt, 1, function(x) sum(sign(x))/length(x))
  avg_abun <- apply(dt, 1, mean)
  dt2 <- data.frame(dt, freq = freq, avg_abun = avg_abun)
  data <- dt2[dt2$freq > low_freq & dt2$avg_abun > low_abun, ]
  data$freq <- NULL
  data$avg_abun <- NULL
  return(data)
}

na.prop <- function(vec) sum(is.na(vec)) / length(vec)

cv <- function(vec) sd(vec)/mean(vec)


na.replace <- function(vec){
  #fill with means
  lev <- levels(factor(vec))
  if(length(lev) > 2){
    vec[is.na(vec)] <- mean(vec, na.rm = T)
  }else{
    #fill with frequency
    vec[is.na(vec)] <- frequency(vec)
  }
  return(vec)
}

###find the index of binary var in a mat
find_binary_var_index <- function(mat, levels=2){
  n.vars <- ncol(mat)
  res <- rep(NA, n.vars)
  for(i in 1:n.vars){
    if(length(levels(factor(mat[, i]))) <= levels){
      res[i] <- TRUE
    }else{
      res[i] <- FALSE
    }
  }
  return(res)
}

###choose top taxa by mean or median or some else
top_choose <- function(mat, margin = 1, FUN = "mean", top = 10){
  mat <- as.data.frame(mat)
  if(nrow(mat) < top)
    stop("The row number of input data must be larger than ", top, "\n")
  mat$tmp <- apply(mat, margin, FUN)
  mat.ord <- mat[order(mat$tmp, decreasing = TRUE), ][1:top,]
  mat.ord$tmp <- NULL
  return(mat.ord)
}

###
w.res.get <- function(var, restype = c('statistics', 'p.value'), 
                      binary_mat, stat_mean = FALSE, adjust = TRUE){
  restype <- match.arg(restype)
  n.bm <- ncol(binary_mat)
  w <- p <- diff.type <- diff.res <- vector(length = n.bm)
  # vars <- rep(deparse(substitute(var)), n.bm)
  for(i in seq_along(1:n.bm)){
    w.res <- wilcox.test(var ~ binary_mat[,i], exact = T)
    w[i] <- w.res$statistic
    p[i] <- w.res$p.value
    mean.res <-  aggregate(var ~ binary_mat[,i], FUN = mean)
    diff.type[i] <- paste0(mean.res[1,1],'-',mean.res[2,1])
    diff.res[i] <- mean.res[1,2] - mean.res[2,2]
  }
  q = p.adjust(p, 'BH')
  diff_R <- data.frame(var = colnames(binary_mat), diff.type, diff.res, w, p, q)
  if(stat_mean){
    return(diff_R)
  }else{
    if(restype == 'statistics'){
      return(w)
    }else if(restype == 'p.value'){
      if(adjust){
	  	return(q)
	  }else{
		return(p)
	  }
    }else{
      return(0)
    }
  }
}

stat_mat_get <- function(micro, binary_mat){
  WMAT <- apply(micro, 2, w.res.get, restype = 'statistics', binary_mat = binary_mat)
  PMAT <- apply(micro, 2, w.res.get, restype = 'p.value', binary_mat = binary_mat)
  rownames(WMAT) <- colnames(binary_mat)
  rownames(PMAT) <- colnames(binary_mat)
  return(list(WMAT = WMAT, PMAT = PMAT))
}

clr <- function(vec){
  # in a sample
  GM <- exp(mean(log(vec), na.rm = T))
  return(log( vec / GM, base = 10))
}

#
getRandomColor <- function(n){
  ###mannual color
  seedcolor <- c('darkorange', 'black', 'darkmagenta', 'darkgreen','cyan4', 'springgreen', 'indianred1', 'red3', 
                 'dodgerblue', 'violetred1', 'turquoise1', 'tan1', 'olivedrab2', 'blue', 'darkgreen', 'red')
  ### random colors()
  all <- grDevices::colors()
  randomcolor <- as.vector(all[grep('grey|gray|white', all, invert = T)])
  # cat(randomcolor)
  set.seed(10)
  if(n <= length(seedcolor)){
    col <- sample(seedcolor, n)
  }else{
    col <- sample(randomcolor, n)
  }
  return(col)
}
