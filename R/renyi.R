
#' Generalized Renyi Transform
#'
#' A Generalization of Aldous Renyi's representation of exponential order statistics
#'
#' Maps a vector of shifted and scaled independent exponential random variables to a sequence of standard independent exponential random variables based on the gaps (jumps) between the initial random variables
#'
#' @references
#' Christ, R., Hall, I. and Steinsaltz, D.  (2024) "The Renyi Outlier Test", \href{https://arxiv.org/abs/2411.13542}{arXiv:2411.13542} . Available at: \doi{10.48550/arXiv.2411.13542}.
#'
#' @param x a vector of independent exponential random variables of the form \eqn{X_j = \eta_j Y_j + \zeta_j} where each \eqn{X_j} is an independent exponential random variable with rate 1
#' @param eta vector of scale parameters implicit in the construction of \code{x}: \code{eta[j]} = \eqn{\eta_j}
#' @param zeta vector of shift parameters implicit in the construction of \code{x}: \code{zeta[j]} = \eqn{\zeta_j}
#' @return a list containing two elements
#' \describe{
#'   \item{`exps`}{a vector of independent standard exponentials where \code{exps[1]} is the exponential jump corresponding to \code{min(x)} and \code{tail(exps,1)} is the exponential jump corresponding to \code{max(x)}.}
#'   \item{`order`}{\code{order(x)}.}
#' }
#' @examples
#' # example code
#'
#' a <- rchisq(10,1)
#' b <- rnorm(10)
#' xx <- a*rexp(10)+b
#' generalized_renyi_transform(xx, a, b)
#' @export
generalized_renyi_transform <- function(x, eta = NULL, zeta = NULL){
  # maps independent exponentials of the form \eta * Y + \zeta where each Y_j is a standard exponential
  # to iid exponentials
  # for the result returned the first element corresponds to exponential jump of the min \eta * Y + \zeta
  # the last element corresponds to the exponential jump of the max \eta * Y + \zeta

  # if a is NULL it is treated as a vector of 1s

  if(is.null(eta) & is.null(zeta)){
    # do standard renyi transformation and exit
    ox <- order(x)
    x <- x[ox]
    return(list("exps" = c(x[1],diff(x)) * seq.int(length(x),1),
                "order" = ox))
  }

  if(is.null(zeta)){
    zeta <- double(length(x))
  }

  p <- length(x)
  if(p == 1){return(list("exps" = if(is.null(eta)){x-zeta}else{(x-zeta)/eta},
                         "order" = 1L))}

  # the algorithm below assumes that min(zeta) = 0, so we shift all of the exponentials
  # so that the first baseline (knot) is at zero
  min_b <- min(zeta)
  x <- x - min_b
  zeta <- zeta - min_b

  if(!is.null(eta)){
    eta_inv <- 1/eta
  }

  ox <- order(x)
  ob <- order(zeta)

  std_expos <- double(p)

  z <- x[ox[1]] # current distance FROM current_baseline (current_baseline initialized at 0)
  running_log_sum <- 0

  j <- 0 # current number of order statistics x we've already PASSED (currently working on the j+1th order statistic)
  k <- 0 # current number of baselines we've already PASSED

  current_expo_ind <- logical(p) # start from none
  next_baseline <- zeta[ob[1]] # FROM current_baseline (current_baseline initialized at 0), so next_baseline is initialized at min(zeta) = the first potential point

  sum_current_expo_ind <- 0
  sum_current_eta_inv <- 0

  while(is.finite(z)){
    if(z > next_baseline){
      z <- z - next_baseline

      if(is.null(eta)){
        running_log_sum <- running_log_sum + next_baseline * sum_current_expo_ind
      } else {
        running_log_sum <- running_log_sum + next_baseline * sum_current_eta_inv
      }

      k <- k+1

      #current_expo_ind[ob[k]] <- TRUE

      if(!current_expo_ind[ob[k]]){
        current_expo_ind[ob[k]] <- TRUE
        sum_current_expo_ind <- sum_current_expo_ind + 1
        if(!is.null(eta)){
          sum_current_eta_inv <- sum_current_eta_inv + eta_inv[ob[k]]
        }
      }

      next_baseline <- if(k == p){
        Inf
      }else{
        zeta[ob[k+1]] - zeta[ob[k]]
      }

    } else {
      # store std_expo
      if(is.null(eta)){
        std_expos[j+1L] <- running_log_sum + z * sum_current_expo_ind # desired output (independent standard exponentials)
      } else {
        std_expos[j+1L] <- running_log_sum + z * sum_current_eta_inv # desired output (independent standard exponentials)
      }
      running_log_sum <- 0

      next_baseline <- next_baseline - z

      j <- j+1

      #current_expo_ind[ox[j]] <- FALSE

      if(current_expo_ind[ox[j]]){
        current_expo_ind[ox[j]] <- FALSE
        sum_current_expo_ind <- sum_current_expo_ind - 1
        if(!is.null(eta)){
          sum_current_eta_inv <- sum_current_eta_inv - eta_inv[ox[j]]
        }
      }


      z <- if(j == p){
        Inf
      } else {
        x[ox[j+1]] - x[ox[j]]
      }
    }

    # if(z==0){
    #   # find next representable double just above x[ox[j]] and divide by 2
    #   # and then simulate a uniform as a difference between 0 and that
    #   x[ox[j]]
    # }


    #if(j == (p-1)){browser()}

  }

  list("exps" = std_expos,
       "order" = ox)
}


# the generalized_renyi_transform reduces to the classic
# renyi transform under eta=1,zeta=0
################################################################
# p <- 1e3
# x <- rexp(p)
# ex <- sort(x)
# all.equal(generalized_renyi_transform(x)$exps,
#           c(ex[1],diff(ex)) * seq.int(length(ex),1))

# test the generalized_renyi_transform against a simple case
# where it's easy to calculate the solution by hand
################################################################

# all.equal(generalized_renyi_transform(
#   x = c(6.5,6.8,7.2,8),
#   zeta = c(4,5,6,7))$exps,
#   c(1+2+3*0.5,2*0.3,1*0.2+2*0.2,1*0.8))
#
# all.equal(generalized_renyi_transform(
#   x = c(6.5,6.8,7.2,8),
#   zeta = c(4,5,5,7))$exps,
#   c(1+3*1.5,2*0.3,1*0.2+2*0.2,1*0.8))
#
# all.equal(generalized_renyi_transform(
#   x = c(6.2,6.4,6.6,8),
#   zeta = c(4,5,6,7))$exps,
#   c(1+2*1+3*0.2,2*0.2,1*0.2,1*1))
#
# the generalized_renyi_transform handles
# gaps where there are no more competing expoentials until
# we hit the next baseline well.
################################################################
# generalized_renyi_transform(c(2,10),eta = c(1,1),zeta = c(0,5))


#' Renyi Outlier Test
#'
#' A fast, numerically precise outlier test for a vector of exact p-values allowing for prior information
#'
#' The about which p-values are outlying and "how much" of an outlier they are expected to be
#'
#' @references
#' Christ, R., Hall, I. and Steinsaltz, D.  (2024) "The Renyi Outlier Test", \href{https://arxiv.org/abs/2411.13542}{arXiv:2411.13542} . Available at: \doi{10.48550/arXiv.2411.13542}.
#'
#' @param u a vector of p-values
#' @param k a rough upper bound on the number of outliers expected to be present in u
#' @param pi optional vector such that \code{pi[j]} is proportional to the probability that \code{u[j]} is an outlier. The default, \code{NULL}, corresponds to \code{pi = rep_len(1,length(u))}.
#' @param eta optional vector proportional to how far outlying we expect \code{u[j]} to be given \code{u[j]} is an outlier.
#' More precisely, in the common context where each element of u can be thought of as a p-value for testing whether some coefficient \eqn{\beta} in a linear regression model is zero, we assume \code{eta[j]} is proportional to \eqn{\mathbb{E}\left[\left. \beta_j^2 \right| \beta_j \neq 0\right]}.
#' The default, \code{NULL}, corresponds to \code{eta = rep_len(1,length(u))}.
#' @return a list containing three elements
#' \describe{
#'   \item{`p_value`}{the p-value returned by the Renyi Outlier Test;}
#'   \item{`max_k`}{a power of 2 in 2^(0:\code{k}) denoting the number of tail p-values that yielded the most significant signal when running the Renyi Outlier Test;}
#'   \item{`p_value_k1`}{the p-value that would be returned by the Renyi Outlier Test assuming \code{k}=1;}
#'   \item{`exit_status`}{a character string describing any problems that may have been encountered during evaluation, "default is no problems";}
#'   \item{`u`}{the vector of p-values used by the outlier test after adjusting the \code{u} provided for \code{pi} and \code{eta}.}
#' }
#' @examples
#' # example code
#'
#' p <- 1e4
#' u <- runif(p)
#' u[c(53,88,32)] <- 1e-6 # add a few outliers
#' renyi(u)$p_value # test for outliers without any prior knowledge
#' renyi(u,pi=c(rep(1,100),rep(10^-3,p-100)))$p_value # test for outliers with prior knowledge
#' @export
renyi <- function(u, k = ceiling(0.01*length(u)), pi = NULL, eta = NULL){

  if(!is.vector(u) || any(!is.finite(u)) || !all(u>=0) || !all(u <= 1)){
    stop("u must be a vector, assumed to be independent standard uniform r.v.s, with all entries in [0,1].")
  }

  if(match(0, u, nomatch = 0L)){
    warning("u passed to renyi contained at least one 0")
    return(list("p_value" = 0,
                "exit_status" = "zero detected in input u, setting output u to all NAs",
                "u" = rep(NA_real_,length(u))))
  }

  p <- length(u)

  if(!is.vector(k) || length(k)!=1 || !is.finite(k) || k < 1 || k!=as.integer(k)){
    stop("k must be a positive integer.")
  }

  k <- 2^(max(0,ceiling(log2(min(128,k,p)))))

  if(!is.null(pi) && (!is.vector(pi) || length(pi) != p || any(!is.finite(pi)) || !all(pi>0) )){
    stop("when provided, pi must be a vector of strictly positive weights proportional to prior probabilities.")
  }

  if(!is.null(eta) && (!is.vector(eta) || length(eta) != p || any(!is.finite(eta)) || !all(eta>0) )){
    stop("when provided, eta must be a vector of strictly positive weights")
  }


  eta_x <- if(is.null(eta)){-log(u)}else{-eta * log(u)}

  if(is.null(pi)){ # we take pi = 1
    rt <- generalized_renyi_transform(eta_x, eta = eta)
  } else {
    zeta <- if(is.null(eta)){log(pi)}else{eta * log(pi)}
    rt <- generalized_renyi_transform(eta_x + zeta, eta = eta, zeta = zeta)
  }

  new_u <- exp(-cumsum(rt$exps / seq.int(p,1))) # sorted from largest u to smallest u

  uk_exp <- -stats::pbeta(new_u[p-k+1L],k,p-k+1L,lower.tail = TRUE,log.p = TRUE)

  if(k==1){
    p_value <- exp(-uk_exp)
    max_k <- 1
    p_value_k1  <- p_value
  } else {
    p_value <- mpse_test(c(rev(utils::tail(rt$exps,k-1L)),uk_exp))
    max_k <- attr(p_value,"max_k")
    attributes(p_value) <- NULL
    p_value_k1 <- exp(stats::pbeta(new_u[p],1,p,lower.tail = TRUE,log.p = TRUE))
  }

  new_u[rt$order] <- new_u

  list("p_value" = p_value,
       "max_k" = max_k,
       "p_value_k1" = p_value_k1,
       "exit_status" = "no problems",
       "u" = new_u) # the effective independent standard uniforms used in calculating the p-value after adjusting for pi and eta
}


# Demonstrate list elements returned by renyi
################################################################
# renyi(runif(1e3),k = 32)


# rd_res <-replicate(5e3,renyi:::mpse_test(rexp(16)))
# shapiro.test(qnorm(rd_res))
# t.test(qnorm(rd_res))


# the renyi returns u provided
# under defaults eta=1,zeta=0.
################################################################

# p <- 1e3
# u <- runif(p)
# res <- renyi(u)
# all.equal(u,res$u)


# the renyi returns uniforms that are very close to u provided
# as we approach the defaults eta=1,zeta=0.
################################################################

# p <- 1e3
# u <- runif(p)
# res <- renyi(u)
# all.equal(u,res$u)
#
# res <- renyi(c(u,runif(p)),pi = c(rep(1,2*p)))
# all.equal(u,res$u[1:p])
#
# res <- renyi(c(u,runif(p)),pi = c(rep(1,p),rep(1e-3,p)))
# plot(u,2*res$u[1:p])
# abline(0,1)


# the renyi matches rdistill::renyi_test
# under defaults eta=1,zeta=0.
################################################################

# p <- 1e3
# u <- runif(p)
# res <- renyi(u)
# all.equal(res$p_value,rdistill::renyi_test(u,k=ceiling(0.01*length(u)))$p_value)


# the renyi is calibrated
################################################################
#
# p <- 1e3
# pi <- runif(p)
#
# n_reps <- 5e3
# res <- data.frame(test_pval = double(n_reps),
#                   st_pval = double(n_reps),
#                   tt_pval = double(n_reps))
#
# for(i in 1:n_reps){
#   u <- runif(p)
#   temp_res <- renyi(u, pi = pi)
#   res$test_pval[i] <- temp_res$p_value
#   zz <- qnorm(temp_res$u)
#   res$st_pval[i] <- shapiro.test(zz)$p.value
#   res$tt_pval[i] <- t.test(zz)$p.value
#   print(i)
# }
#
# plot(u,temp_res$u)
#
# apply(res,2,function(x){shapiro.test(qnorm(x))$p.value})
# apply(res,2,function(x){t.test(qnorm(x))$p.value})
#
# # st_pvalues are not exact p-values -- slightly conservative (subuniform) -- which is fine and why they fail the global test
# qqnorm(qnorm(res$st_pval))
# qqline(qnorm(res$st_pval))
# #



# Demonstrate how that renyi
# adjusts the multiple testing burden
################################################################

# p <- 1e4
#
# u <- runif(p)
# all.equal(u,renyi(u)$u)
#
# u <- data.frame("helpful_pi" = runif(p),
#                 "unhelpful_pi" = runif(p))
#
# u$helpful_pi[c(53,88,32)] <- 1e-6
# u$unhelpful_pi[c(530,880,320)] <- 1e-6
#
#
#
# xx <- 0:4
# res <- as.data.frame(matrix(0,length(xx),ncol(u)))
#
# for(j in 1:ncol(u)){
#   for(i in 1:length(xx)){
#     res[i,j] <- -log10(renyi(u[,j], pi = c(rep(1,100),rep(10^-xx[i],p-100)), k = 4)$p_value)
#   }
# }
#
# plot(xx, res[,1], ylim=c(0,12),col="darkorange",type="l",las=1,bty="n",ylab="-log10 p-value",xlab="-log10 pi",lwd=3)
# abline(h=-log10(renyi(u$helpful_pi[1:100], k = 4)$p_value),lty=2,lwd=3)
# lines(xx, res[,2],col="skyblue",lwd=3)








