#' Exact hazard rate test - single sample
#' 
#' Performs an exact test of the null that a single sample hazard rate is equal to some value.
#' 
#' @param time a numeric vector of survival times.
#' @param status status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param null the null hazard rate. Default is 1.
#' @param conf.level confidence level of the returned confidence interval. Must be a single number between 0 and 1.
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'    \item{p.value}{The p-value of the test.} 
#'    \item{estimate}{An estimate of the hazard rate.} 
#'    \item{null.value}{The null hazard rate.} 
#'    \item{conf.int}{A confidence interval for the hazard rate.} 
#'    \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided".} 
#'    \item{method}{A character string indicating the method employed.} 
#'    \item{data.name}{A character string giving the name of the data.} 
#'   }
#' @references Newman (2001), page 203-4.
#' @examples 
#' ## Example 10.3
#' exact.rate.test(time = 10, status = 2, null = 0.4)
#' 
#' ## With Breast Cancer Survival data
#' with(breast.survival, exact.rate.test(time = time, status = status))
#' with(breast.survival, exact.rate.test(time = time, status = status, null = 0.01))
exact.rate.test <- function(time, status, null=1, conf.level = 0.95){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  # estimate
  d <- sum(status, na.rm = TRUE)
  n <- sum(time, na.rm = TRUE)
  est <- d/n
  names(est) <- "hazard rate"
  names(null) <- names(est)
  # CI
  alpha <- (1-conf.level)/2
  f1 <- function(x)1 - ppois(q = d-1, lambda = n*x) - alpha
  f2 <- function(x)ppois(q = d, lambda = n*x) - alpha
  f1.out <- uniroot(f1, interval = c(0,n))
  f2.out <- uniroot(f2, interval = c(0,n))
  CINT <- c(f1.out$root, f2.out$root)
  attr(CINT, "conf.level") <- conf.level
  
  # test
  p.value <- min(ppois(q = d, lambda = n*null) * 2, 1)
  RVAL <- list(p.value = p.value,estimate = est, null.value = null,
               conf.int = CINT, alternative = alternative,
               method = "Exact Hazard Rate Test for a single sample", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

