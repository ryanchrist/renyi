% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/renyi.R
\name{renyi}
\alias{renyi}
\title{Renyi Outlier Test}
\usage{
renyi(u, k = ceiling(0.01 * length(u)), pi = NULL, eta = NULL)
}
\arguments{
\item{u}{a vector of p-values}

\item{k}{a rough upper bound on the number of outliers expected to be present in u}

\item{pi}{optional vector such that \code{pi[j]} is proportional to the probability that \code{u[j]} is an outlier. The default, \code{NULL}, corresponds to \code{pi = rep_len(1,length(u))}.}

\item{eta}{optional vector proportional to how far outlying we expect \code{u[j]} to be given \code{u[j]} is an outlier.
More precisely, in the common context where each element of u can be thought of as a p-value for testing whether some coefficient \eqn{\beta} in a linear regression model is zero, we assume \code{eta[j]} is proportional to \eqn{\mathbb{E}\left[\left. \beta_j^2 \right| \beta_j \neq 0\right]}.
The default, \code{NULL}, corresponds to \code{eta = rep_len(1,length(u))}.}
}
\value{
a list containing three elements
\describe{
  \item{`p_value`}{the p-value returned by the Renyi Outlier Test;}
  \item{`exit_status`}{a character string describing any problems that may have been encountered during evaluation, "default is no problems";}
  \item{`u`}{the vector of p-values used by the outlier test after adjusting the \code{u} provided for \code{pi} and \code{eta}.}
}
}
\description{
A fast, numerically precise outlier test for a vector of exact p-values allowing for prior information
}
\details{
The about which p-values are outlying and "how much" of an outlier they are expected to be
}
\examples{
# example code

p <- 1e4
u <- runif(p)
u[c(53,88,32)] <- 1e-6 # add a few outliers
renyi(u)$p_value # test for outliers without any prior knowledge
renyi(u,pi=c(rep(1,100),rep(10^-3,p-100)))$p_value # test for outliers with prior knowledge
}
\references{
Christ, R., Hall, I. and Steinsaltz, D.  (2024) "The Renyi Outlier Test", \href{https://arxiv.org/abs/2411.13542}{arXiv:2411.13542} . Available at: \doi{10.48550/arXiv.2411.13542}.
}
