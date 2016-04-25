#' Test for linear trend (in log-hazard ratios)
#' 
#' Performs a test for linear trend of log-hazard ratios in stratified survival data.
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param exposure a factor vector with two levels indicating exposure. The 
#'   first level is assumed to be the exposed condition.
#' @param strata a factor vector with at least two levels indicating strata.
#' @param score group score.
#' @return An object of class "htest" with title, test statistic, and p-value.
#' @references Newman (2001), page 221.
#' @examples 
#' ## Example 10.14
#' with(breast.survival, linear.trend.test(time, status, receptor.level, stage)) 
#' 
#' ## Using specified scores
#' with(breast.survival, linear.trend.test(time, status, receptor.level, stage,
#'                                         scores = c(1,2,6))) 

linear.trend.test <- function(time, status, exposure, strata, 
                              scores = seq_along(length(levels(strata)))){
  dname <- paste("\n  ", deparse(substitute(time)), 
                 "\n  ", deparse(substitute(status)),
                 "\n   using scores:", paste(scores, collapse = " "))
  alternative <- "two.sided"
  
  dk <- tapply(status, list(exposure, strata), sum)
  nk <- tapply(time, list(exposure, strata), sum)
  mk <- apply(dk,2,sum)
  n <- apply(nk,2,sum)
  # estimate common hazard ratio; need this for test statistic calculation
  f1 <- function(x){
    sum((x*mk*nk[1,])/(x*nk[1,] + nk[2,])) - sum(dk[1,])
  }
  est <- uniroot(f = f1, interval = c(0,1e5))$root
  hr2 <- mk/(est*nk[1,] + nk[2,])
  hr1 <- est*hr2
  # fitted counts:
  dh1 <- hr1*nk[1,]
  dh2 <- hr2*nk[2,]
  s <- scores
  num <- sum(s * (dk[1,] - dh1))^2
  vk <- (1/dh1 + 1/dh2)^(-1)
  den <- sum(s^2 * vk) - (sum(s * vk))^2 / sum(vk)
  STATISTIC <- num/den
  p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  names(STATISTIC) <- "X-squared"
  METHOD <- paste("Test for linear trend (in log-hazard ratios)")
  RVAL <- list(statistic = STATISTIC, parameter = c(df = 1), p.value = p.value,
               method = METHOD, 
               data.name = dname)  
  class(RVAL) <- "htest"
  return(RVAL) 
  
}
