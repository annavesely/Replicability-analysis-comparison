
require(hommel)
require(adaFilter)


# -------------------------------------------------------


# Unconditional PC p-value through Fisher combination
# x: vector of sorted p-values
# gamma: replicability level
# output: PC p-value

poolPvals <- function(x, gamma){
  z <- x[gamma : length(x)]
  comb <- -2*sum(log(z))
  out <- pchisq(comb, df=length(z)*2, lower.tail=FALSE)
  return(out)
}


# -------------------------------------------------------


# Benjamini, Heller, Yekutieli (2009)
# pmat0: matrix of pvalues (m rows = features, s columns = studies);
# alpha: significance level
# output: boolean vector of length m, keeping track of replicability

BHY <- function(pmat0, alpha=0.05){
  
  pmat <- t(apply(pmat0, 1, function(x) sort(x)))
  m <- nrow(pmat)
  s <- ncol(pmat)
  
  X <- rep(0,m)
  
  pv <- apply(pmat, 1, function(x) poolPvals(x, 1))
  tmp <- p.adjust(pv, method="BH")
  sel <- which(tmp <= alpha)
  m1 <- length(sel)
  thr <- alpha * m1 / m
  
  for(j in sel){
    gamma <- 0
    
    while(gamma < s){
      
      gamma <- gamma + 1
      
      pv <- poolPvals(pmat[j,], gamma)
      if(pv > thr){
        gamma <- gamma - 1
        break
      }
    }
    X[j] <- gamma
  }
  
  return(X)
}


# -------------------------------------------------------


# adaFilter: Wang, Gui, Su, Sabatti, Owen (2020)
# pmat0: matrix of pvalues (m rows = features, s columns = studies);
# alpha: significance level
# output: boolean vector of length m, keeping track of replicability

adafilter <- function(pmat0, alpha=0.05){
  
  m <- nrow(pmat0)
  s <- ncol(pmat0)
  
  X <- rep(0,m)
  
  for(gamma in (2:s)){
    rej <- adaFilter(pmat0, r=gamma, type.I.err="FDR", alpha=alpha)$decision
    rej <- as.logical(rej)
    X[rej] <- gamma
  }
  
  return(X)
}


# -------------------------------------------------------


# ARI: Goeman and Solari (2011)
# pmat0: matrix of pvalues (m rows = features, s columns = studies);
# alpha: significance level
# output: boolean vector of length m, keeping track of replicability

ari <- function(pmat0, alpha=0.05){
  
  m <- nrow(pmat0)
  s <- ncol(pmat0)
  
  p <- as.vector(pmat0)
  hom <- hommel::hommel(p)
  
  out <- rep(0,m)
  
  for(j in 1:m){
    ix <- (0:(s-1))*m + j
    out[j] <- hommel::discoveries(hom, ix=ix, alpha=alpha)
  }
  
  return(out)
}


