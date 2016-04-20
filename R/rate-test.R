#' Asymptotic hazard rate test - single sample
#' 
#' Performs an asymtotic test of the null that a single sample hazard rate is equal to some value.
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param null the null hazard rate. Default is 1. 
#' @param conf.level confidence level of the returned confidence interval. Must be a single number between 0 and 1.
#' @param explicit a logical indicating whether to calculate an explicit or implicit confidence interval. Default is TRUE.
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'    \item{statistic}{The chi-square test statistic.} 
#'    \item{p.value}{The p-value of the test.} 
#'    \item{estimate}{An estimate of the hazard rate.} 
#'    \item{null.value}{The null hazard rate.} 
#'    \item{conf.int}{A confidence interval for the hazard rate.} 
#'    \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided".} 
#'    \item{method}{A character string indicating the method employed.} 
#'    \item{data.name}{A character string giving the name of the data.} 
#'   }
#' @references Newman (2001), page 205.
#' @examples 
#' ## Examples 10.5 and 10.6
#' rate.test(time = 10, status = 2, null = 0.4)
#' rate.test(time = 25, status = 5, null = 0.4)
#' rate.test(time = 50, status = 10, null = 0.4)
#' 
#' rate.test(time = 10, status = 2, null = 0.4, explicit = FALSE)
#' rate.test(time = 25, status = 5, null = 0.4, explicit = FALSE)
#' rate.test(time = 50, status = 10, null = 0.4, explicit = FALSE)
#' 
#' ## Example 10.7
#' with(breast.survival, rate.test(time = time, status = status, explicit = FALSE))
#' 
rate.test <- function(time, status, null=1, conf.level = 0.95, explicit=TRUE){
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
  if(explicit){
    CINT <- est + c(-1,1)*((qnorm(1 - alpha)*sqrt(d))/n)  
  } else {
    a <- n^2
    b <- -n*(2*d + qnorm(1 - alpha)^2)
    c <- d^2
    CINT <- (-b + c(-1,1)*sqrt(b^2 - 4*a*c))/(2*a)
  }
  attr(CINT, "conf.level") <- conf.level
  
  # test
  STATISTIC <- ((d - null*n)^2)/(null*n)
  names(STATISTIC) <- "X-squared"
  p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, 
               null.value = null,
               conf.int = CINT, alternative = alternative,
               method = "Asymptotic Hazard Rate Test for a single sample", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

