#' Likelihood ratio test of homogeneity of hazard ratios (stratified survival 
#' data)
#' 
#' Performs a likelihood ratio test of homogeneity of hazard ratios for 
#' stratified survival data. The null hypothesis is that all stratum-specific 
#' hazard ratios are equal.
#' 
#' According to Newman (p. 124): "the homogeneity assumption is merely a 
#' convenient fiction that is adopted in order to simplify the analysis and 
#' interpretation of data." A suggested approach is to first examine 
#' stratum-specific hazard ratios and their confidence intervals to get a sense 
#' if there are meaningful differences across strata (see 
#' \code{\link{hazard.ratios}}); perform a formal test of homogeneity; and then 
#' "synthesize this information along with substantive knowledge of the 
#' relationship between exposure and disease, taking into account the aims of
#' the study."
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param exposure a factor vector with two levels indicating exposure. The 
#'   first level is assumed to be the exposed condition.
#' @param strata a factor vector with at least two levels indicating strata.
#' @return An object of class "htest" with title, test statistic, and p-value.
#' @export
#' @references Newman (2001), page 221.
#' @examples 
#' ## Example 10.14
#' with(breast.survival, lrt.homogeneity(time, status, receptor.level, stage)) 
lrt.homogeneity <- function(time, status, exposure, strata){
  dname <- paste("\n  ", deparse(substitute(time)), 
                 "\n  ", deparse(substitute(status)))
  dk <- tapply(status, list(exposure, strata), sum)
  nk <- tapply(time, list(exposure, strata), sum)
  mk <- apply(dk,2,sum)
  n <- apply(nk,2,sum)
  # estimate common hazard ratio; need this for test statistic calculation
  f1 <- function(x){
    sum((x*mk*nk[1,])/(x*nk[1,] + nk[2,])) - sum(dk[1,])
  }
  est <- stats::uniroot(f = f1, interval = c(0,1e5))$root
  hr2 <- mk/(est*nk[1,] + nk[2,])
  hr1 <- est*hr2
  # fitted counts:
  dh1 <- hr1*nk[1,]
  dh2 <- hr2*nk[2,]
  
  STATISTIC <- 2 * sum(dk[1,] * log((dk[1,]/dh1)) + dk[2,] * log((dk[2,]/dh2)))
  df <- ncol(dk) - 1
  names(df) <- "df"
  p.value <- stats::pchisq(STATISTIC, df = df, lower.tail = FALSE)
  names(STATISTIC) <- "X-squared"
  METHOD <- paste("Likelihood ratio test of homogeneity")
  RVAL <- list(statistic = STATISTIC, parameter = df, p.value = p.value,
               method = METHOD, 
               data.name = dname)  
  class(RVAL) <- "htest"
  return(RVAL) 
  
}

