#' Hazard ratios of stratified survival data
#' 
#' Returns hazard ratios and confidence intervals for each stratum.
#' 
#' Hazard ratios and their confidence intervals are calculated with the Wald 
#' method using \code{epitools::rateratio}. This can be changed to 
#' median-unbiased estimation by setting \code{method="midp"}.
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param exposure a factor vector with two levels indicating exposure. The 
#'   first level is assumed to be the exposed condition.
#' @param strata a factor vector with at least two levels indicating strata.
#' @param method method for calculating hazard ratio and confidence interval.
#'   Choices are "wald" (default) and "midp".
#' @return  A list containing a hazard ratio estimate and confidence interval for each stratum.
#' @export
#' @references Newman (2001), page 208.
#' @examples 
#' ## Table 10.15 from Example 10.14
#' with(breast.survival, hazard.ratios(time, status, receptor.level, stage))
#' 
#' ## using the "midp" method
#' with(breast.survival, hazard.ratios(time, status, receptor.level, stage, method="midp"))
hazard.ratios <- function(time, status, exposure, strata, method="wald"){
  dk <- tapply(status, list(exposure, strata), sum)
  nk <- tapply(time, list(exposure, strata), sum)
  est <- lapply(seq(ncol(dk)), 
                function(x,y=method)epitools::rateratio(x = dk[,x],y = nk[,x],
                                               method = y, 
                                               rev = "columns")$measure[2,])
  names(est) <- levels(strata)
  return(est)
}
