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
#' @param group a factor vector with at least two levels
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
#' @references Newman (2001), page 216.
#' @examples 
#' with(breast.survival, mantelhaen.pt.test(time, status, stage))

mantelhaen.pt.test <- function(time, status, group){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  null <- 1
  d <- tapply(status, group, sum)
  n <- tapply(time, group, sum)
  est <- (d*n[1])/(d[1]*n)
  names(est) <- levels(group)
  names(null) <- "hazard ratio"
  e <- n * sum(d) / sum(n)
  STATISTIC <- sum((d - e)^2 / e)
  p.value <- pchisq(STATISTIC, df = length(d) - 1, lower.tail = FALSE)
  if(any(e < 5)) warning(paste("expected counts less than 5:",e[e < 5]))
  names(STATISTIC) <- "Mantel-Haenszel X-squared"
  p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, null.value = null,
               alternative = alternative,
               method = "Mantel-Haenszel test of association for person-time data", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}
