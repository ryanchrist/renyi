# copied over from QForm
eval.cdf.pos <- function(q, cdf, cdf.body, density = FALSE, lower.tail = TRUE, log.p = FALSE){

  if(density){

    if(log.p){

      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( log(-cdf$b.l)-(cdf$a.l + (cdf$b.l+1)*log(q)) ),
                                     ifelse(q>cdf$x[cdf$n],log(cdf$b.r)-(cdf$a.r + cdf$b.r*q),
                                            log(cdf.body(q,deriv=1))))))
    }else{

      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( -cdf$b.l*exp(-(cdf$a.l + (cdf$b.l+1)*log(q))) ),
                                     ifelse(q>cdf$x[cdf$n],cdf$b.r*exp(-(cdf$a.r + cdf$b.r*q)),
                                            cdf.body(q,deriv=1)))))
    }

  }else{

    if(lower.tail & !log.p){
      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( exp(-(cdf$a.l + cdf$b.l*log(q))) ),
                                     ifelse(q>cdf$x[cdf$n],suppressWarnings( -expm1(-(cdf$a.r + cdf$b.r*q)) ),
                                            cdf.body(q)))))
    }

    if(!lower.tail & !log.p){
      return(ifelse(q < 0, 1, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( -expm1(-(cdf$a.l + cdf$b.l*log(q))) ),
                                     ifelse(q>cdf$x[cdf$n],suppressWarnings( exp(-(cdf$a.r + cdf$b.r*q)) ),
                                            1-cdf.body(q)))))
    }

    if(lower.tail & log.p){
      return(ifelse(q < 0, -Inf, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( -(cdf$a.l + cdf$b.l*log(q)) ),
                                        ifelse(q>cdf$x[cdf$n],suppressWarnings( log(-expm1(-(cdf$a.r + cdf$b.r*q))) ),
                                               suppressWarnings(log(cdf.body(q)))))))
    }
    if(!lower.tail & log.p){
      return(ifelse(q < 0, 0, ifelse(q<cdf$x[1] & (!is.null(cdf$b.l)),suppressWarnings( log(-expm1(-(cdf$a.l + cdf$b.l*log(q))) )),
                                     ifelse(q>cdf$x[cdf$n],-(cdf$a.r + cdf$b.r*q),
                                            suppressWarnings(log1p(-cdf.body(q)))))))
    }
  }
}

# Copied over from QForm with the modification that we only ever use eval.cdf.pos since
# in renyi we're only interested in random variables with non-negative support
wrap.QFcdf <- function(cdf){
  # Returns a function that will evaluate the CDF pointwise
  cdf.body <- stats::splinefun(cdf$x,cdf$y,method="mono")
  if(cdf$type != "pos"){
    stop("cdf must have positive support")
  }
  function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) eval.cdf.pos(x, cdf, cdf.body, density = density, lower.tail = lower.tail, log.p = log.p)
}
