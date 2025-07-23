
internal_mpset <- function(cdf.hat,k,a.r,b.r,x,lower.tail,log.p){
  # Confirm algebraic solution to root matches what we get from machine:
  intersection.x <- (a.r + log(k) ) / (1-b.r)
  #uniroot(function(x){max(x-log(k),0) + cdf.hat(x,density = FALSE,lower.tail = FALSE, log.p = TRUE)},interval = c(5,50))

  y <- cdf.hat(x, density = FALSE, lower.tail = lower.tail, log.p = log.p)

  if(!is.finite(intersection.x) || intersection.x <= log(k)){ # the extrapolated CDF never crosses over the Bonferroni line
    # so we can just return the estimate
    return(y)
  }

  if(length(x) <=0){
    warning("x provided has length 0, returning NA")
    res <- NA_real_
  }

  if(length(x) == 1){
    res <- if(!log.p){
      if(lower.tail){
        if(x > intersection.x){1-k*exp(-x)}else{y}
      } else {
        if(x > intersection.x){k*exp(-x)}else{y}
      }
    } else {
      if(lower.tail){
        if(x > intersection.x){log1p(-k*exp(-x))}else{y}
      } else {
        if(x > intersection.x){-x+log(k)}else{y}
      }
    }
  }

  if(length(x) > 1){
    res <- if(!log.p){
      if(lower.tail){
        ifelse(x > intersection.x,1-k*exp(-x),y)
      } else {
        ifelse(x > intersection.x,k*exp(-x),y)
      }
    } else {
      if(lower.tail){
        ifelse(x > intersection.x,log1p(-k*exp(-x)),y)
      } else {
        ifelse(x > intersection.x,-x+log(k),y)
      }
    }
  }

  res
}


mpset_cdf_2 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- wrap.QFcdf(MPSETestData[["2"]])
  k <- MPSETestData[["2"]]$k
  a.r <- MPSETestData[["2"]]$a.r
  b.r <- MPSETestData[["2"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}

mpset_cdf_4 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- wrap.QFcdf(MPSETestData[["4"]])
  k <- MPSETestData[["4"]]$k
  a.r <- MPSETestData[["4"]]$a.r
  b.r <- MPSETestData[["4"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}

mpset_cdf_8 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- wrap.QFcdf(MPSETestData[["8"]])
  k <- MPSETestData[["8"]]$k
  a.r <- MPSETestData[["8"]]$a.r
  b.r <- MPSETestData[["8"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}

mpset_cdf_16 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- wrap.QFcdf(MPSETestData[["16"]])
  k <- MPSETestData[["16"]]$k
  a.r <- MPSETestData[["16"]]$a.r
  b.r <- MPSETestData[["16"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}

mpset_cdf_32 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- wrap.QFcdf(MPSETestData[["32"]])
  k <- MPSETestData[["32"]]$k
  a.r <- MPSETestData[["32"]]$a.r
  b.r <- MPSETestData[["32"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}

mpset_cdf_64 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- wrap.QFcdf(MPSETestData[["64"]])
  k <- MPSETestData[["64"]]$k
  a.r <- MPSETestData[["64"]]$a.r
  b.r <- MPSETestData[["64"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}

mpset_cdf_128 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- wrap.QFcdf(MPSETestData[["128"]])
  k <- MPSETestData[["128"]]$k
  a.r <- MPSETestData[["128"]]$a.r
  b.r <- MPSETestData[["128"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}


# mpse_test is the same as the original mpse.test
mpse_test <- function(x){
  # x here must be iid standard exponentials, a vector of length 1, 2, 4, 8, 16, 32, 64, or 128

  if(!is.vector(x) || ! length(x) %in% 2^c(0:7)){
    stop("for now, x must be a vector of length 1, 2, 4, 8, 16, 32, 64, or 128")}

  if(length(x)==1){return(stats::pgamma(x,1,lower.tail = FALSE))}

  y <- cumsum(x)[2^(0:log2(length(x)))]

  tm <- 0
  max_k <- 1
  for(i in 0:log2(length(x))){
    new_value <- -stats::pgamma(y[i+1],2^i,lower.tail = FALSE,log.p = TRUE)
    max_idx <- which.max(c(tm, new_value))
    if(max_idx==2){
      tm <- new_value
      max_k <- 2^i
    }
  }

  p_value <- switch(as.character(length(x)),
         "2"   = mpset_cdf_2(tm, lower.tail = FALSE),
         "4"   = mpset_cdf_4(tm, lower.tail = FALSE),
         "8"   = mpset_cdf_8(tm, lower.tail = FALSE),
         "16"  = mpset_cdf_16(tm, lower.tail = FALSE),
         "32"  = mpset_cdf_32(tm, lower.tail = FALSE),
         "64"  = mpset_cdf_64(tm, lower.tail = FALSE),
         "128" = mpset_cdf_128(tm, lower.tail = FALSE))

  attr(p_value,"max_k") <- max_k

  p_value
}



# mks <- replicate(5e3,attr(renyi:::mpse_test(rexp(128)),"max_k"))
#
# barplot(table(mks))
#
#
# ans <- replicate(5e3,renyi:::mpse_test(rexp(4)))
#
# qqnorm(qnorm(ans));abline(0,1)
# shapiro.test(qnorm(ans))
#
#
# ans <- rep(0,1e4)
# for(i in 1:length(ans)){
#   ans[i] <- mpse.test(rexp(64))
#   print(i)
# }
# F.hat <- ecdf(-log10(ans))
# tail.est <- function(x){-log10(1-F.hat(x))}
# xx <- seq(0,5,len=1e3)
# plot(xx,tail.est(xx),type="l")
# abline(0,1,lty=2)
#
# xx <- seq(0,50,len=1e3)
# plot(xx,-log10(rdistill:::mpset_cdf_128(xx,lower.tail = FALSE)),type="l")
