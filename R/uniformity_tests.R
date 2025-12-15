#' Test for Uniformity on [0,1] Using Multiple Statistical Tests
#'
#' A wraper function that performs multiple statistical tests to assess whether a numeric vector represents
#' independent draws from a uniform distribution on the interval [0,1]. The function
#' combines several complementary approaches including tests based on the Rényi Outlier Test (see \code{\link[renyi]{renyi}}),
#' distribution fitting (Kolmogorov-Smirnov), location (t-test), and normality after
#' transformation (Shapiro-Wilk).
#'
#' @param u A numeric vector of values assumed to be on [0,1]. Each element should
#'   represent an independent draw that is being tested for uniformity.
#' @param ... optional arguments to be passed to the Rényi Outlier Test, see \code{\link[renyi]{renyi}}.
#' @return A named list containing p-values from four uniformity tests:
#' \describe{
#'   \item{shapiro}{P-value from Shapiro-Wilk test applied to normal quantile
#'     transformed data. Tests whether Φ⁻¹(u) follows standard normal distribution.}
#'   \item{renyi}{P-value from the Rényi Outlier Test.}
#'   \item{t}{P-value from one-sample t-test testing if mean of Φ⁻¹(u) equals 0.}
#'   \item{ks}{P-value from Kolmogorov-Smirnov test comparing the empirical distribution
#'     to the uniform distribution.}
#' }
#'
#' @details
#' The function applies four different statistical tests:
#'
#' \itemize{
#'   \item \strong{Kolmogorov-Smirnov}: Compares the empirical distribution
#'     of u to the uniform distribution on [0,1]
#'   \item \strong{Rényi Outlier Test}: Tests whether there are outlying small entries of \code{u},see \code{\link[renyi]{renyi}}.
#'   \item \strong{t-test}: Transforms u using the inverse normal CDF and tests
#'     whether the mean equals 0 (expected value for standard normal)
#'   \item \strong{Shapiro-Wilk}: Tests whether the normal quantile transform
#'     Φ⁻¹(u) follows a standard normal distribution
#' }
#'
#' For the Shapiro-Wilk test, if the sample size exceeds 5,000, the function
#' automatically subsamples to 5,000 quantiles to meet the test's sample size
#' limitations.
#'
#' All tests return p-values where small values (typically < 0.05) suggest
#' evidence against the null hypothesis of uniformity.
#'
#' @examples
#' # Test truly uniform data
#' uniform_data <- runif(1000)
#' results1 <- uniformity_tests(uniform_data)
#' print(results1)  # Should show large p-values
#'
#' # Test non-uniform data (beta distribution)
#' beta_data <- rbeta(1000, 2, 5)
#' results2 <- uniformity_tests(beta_data)
#' print(results2)  # Should show small p-values
#'
#' # Test a data with small u outliers
#' outlier_data <- c(uniform_data,1e-5,5e-6,1e-6)
#' results3 <- uniformity_tests(outlier_data)
#' print(results3)  # Should show small p-values
#'
#' # Test while passing a different argument k to the Rényi Outlier Test
#' results4 <- uniformity_tests(outlier_data, k = 4)
#' print(results4)
#'
#'
#' @seealso
#' \code{\link[stats]{ks.test}}, \code{\link[stats]{shapiro.test}},
#' \code{\link[stats]{t.test}}, \code{\link[renyi]{renyi}}
#'
#' @export
uniformity_tests <- function(u, k = min(ceiling(0.05 * length(u)), 32)){

  if (length(u) < 3) {
    stop("'u' must have at least 3 elements.")
  }

  renyi_res <- renyi::renyi(u, k = k)$p_value
  ks_res <- ks.test(u,"punif")$p.value
  z <- qnorm(u)
  t_res <- t.test(z)$p.value

  if(length(u)>5e3){z <- qnorm(quantile(u, probs = ppoints(5e3)))}

  shapiro_res <- shapiro.test(z)$p.value

  list( "shapiro" = shapiro_res,
        "renyi" = renyi_res,
        "t" = t_res,
        "ks" = ks_res)
}
