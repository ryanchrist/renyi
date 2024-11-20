
n <- 1e8
k <- 8
M <- matrix(0,n,k)


M[,1] <- rgamma(n,1)
for(i in 1:7){
  M[,i+1] <- M[,i] + rgamma(n,2^(i-1))
}

M[,1] <- -pgamma(M[,1],1,lower.tail = FALSE,log.p = TRUE)
for(i in 1:7){
  M[,i+1] <- pmax(M[,i],-pgamma(M[,i+1],2^i,lower.tail = FALSE,log.p = TRUE))
  print(paste(i,"th max calculated"))
}


MPSETestData <- as.list(1:7)
names(MPSETestData) <- as.character(2^(1:7))

for(i in 1:7){

  F.true.max <- ecdf(M[,i+1])
  xx <- seq(0,20,len=1e4)
  n <- length(xx)
  fft_cdf <- F.true.max(xx)


  #plot(xx,fft_cdf)

  # Begin search for extrapolation points
  ########################################

  # Focus in on the tail
  ctstart.r <- which.max(fft_cdf>0.99)
  ctstart.l <- which.max(fft_cdf>0.01)-1

  # First: Eliminate Parts of the Estimated CDF that deviate below 0 or above 1 (allowed to hit 0 or 1 at the
  # last point if all of the eigenvalues have the same sign and the domain of the CDF has a bound on that side)

  r0 <- which(fft_cdf[ctstart.r:n] >= 1)

  best.r <- ifelse(length(r0)==0, n, ctstart.r-2+min(r0)) # One point in from the trouble point

  # Second: Eliminate points as unstable if/where the CDF stops being monotonic

  # Take -log density
  log.cdf.r <- suppressWarnings(-log1p(-fft_cdf))

  rblips <- which(diff(log.cdf.r[ctstart.r:best.r]) < -sqrt(.Machine$double.eps))
  best.r <- ifelse(length(rblips)==0,best.r,ctstart.r-2+min(rblips))
  # Again, we take the point one in from the trouble


  # Finally: refine extrapolation points by cropping off the density when it starts becoming curvy due to numerical instability

  limit.r <- Inf
  right.tail <- QForm:::extrapolate_tail(log.cdf.r,xx,ctstart.r,best.r,num.windows=20,right.side=TRUE)

  if(right.tail$successful){
    best.r<-right.tail$best
    b.r<-right.tail$b
    a.r <- log.cdf.r[best.r]-xx[best.r]*b.r
  }

  MPSETestData[[i]] <- list("x" = xx[1:best.r],
                            "y" = fft_cdf[1:best.r],
                            "type"= "pos",
                            "n" = best.r,
                            "interval.width"= NA,
                            "limit.l" = NA,
                            "a.l" = NA,
                            "b.l" = NA,
                            "limit.r" = limit.r,
                            "a.r" = a.r,
                            "b.r" = b.r,
                            "mu" = 0,
                            "Q.sd" = 1,
                            "k" = k)
  print(paste(i,"th dataset calculated"))
}

usethis::use_data(MPSETestData, overwrite = TRUE, internal = TRUE)

xxx <- seq(0,20,len=1e3)
plot(xxx,-log(1-F.true.max(xxx)),col="blue",type="l",
     ylim=c(0,20),xlim=c(0,max(xxx)),bty="n",las=1)
lines(xxx,(a.r + b.r * xxx),col="green")
lines(xxx,pmax(xxx - log(8),0),col="red")
# lines(xxx,
#       -cdf.hat(xxx,lower.tail = FALSE,log.p = TRUE),
#       col="darkorange",lty=2,lwd=2)
abline(0,1,lty=2)
# lines(xxx,
#       -cdf1(xxx,lower.tail = FALSE,log.p = TRUE),
#       col="forestgreen",lty=2,lwd=2)

1/exp(-log(1-F.true.max(8)) - 8)

