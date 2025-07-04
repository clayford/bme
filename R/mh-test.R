#' Mantel-Haenszel Test of association for person-time data
#' 
#' Performs a Mantel-Haenszel test of the null that the hazard ratio is 1 for 
#' all exposure categories, relative to the reference category. A warning is 
#' given if any expected count is less than 5, along with the offending expected
#' count(s).
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param exposure a factor vector with two levels indicating exposure. The
#'   first level is assumed to be the exposed condition.
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'    \item{statistic}{The Mantel-Haenszel X-squared chi-square test statistic.} 
#'    \item{p.value}{The p-value of the test.} 
#'    \item{estimate}{Hazard ratio estimates. The reference level is 1.} 
#'    \item{null.value}{The null hazard ratio, which is currently set to 1.} 
#'    \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided".} 
#'    \item{method}{A character string indicating the method employed.} 
#'    \item{data.name}{A character string giving the name of the data.} 
#'   }
#' @export
#' @references Newman (2001), page 216.
#' @examples 
#' ## Example 10.12
#' with(breast.survival, mantelhaen.pt.test(time, status, stage))
#' 
#' ## Example 10.13, assessment of the Poisson-Exponential Assumption
#' mantelhaen.pt.test(time = c(2363, 7108), status = c(5,44), exposure = gl(n = 2, k = 1))
#' ## p = 0.017; moderate evidence of unequal hazard rates in the two time intervals.

mantelhaen.pt.test <- function(time, status, exposure){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  null <- 1
  d <- tapply(status, exposure, sum)
  n <- tapply(time, exposure, sum)
  est <- (d*n[1])/(d[1]*n)
  names(est) <- levels(exposure)
  names(null) <- "hazard ratio"
  e <- n * sum(d) / sum(n)
  STATISTIC <- sum((d - e)^2 / e)
  p.value <- stats::pchisq(STATISTIC, df = length(d) - 1, lower.tail = FALSE)
  if(any(e < 5)) warning(paste("expected counts less than 5:",e[e < 5]))
  names(STATISTIC) <- "Mantel-Haenszel X-squared"
  p.value <- stats::pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, null.value = null,
               alternative = alternative,
               method = "Mantel-Haenszel test of association for person-time data", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}
