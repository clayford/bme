#' Hazard Ratio Test of association for survival data with k strata
#' 
#' Performs a test of the null that the hazard ratio between \emph{k} strata is 1.
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param exposure a factor vector with at least two levels indicating exposure. The
#'   first level is assumed to be the exposed condition.
#' @param strata a factor vector with at least two levels indicating strata.
#' @param Wald a logical indicating whether to use the Wald or Likelihood Ratio 
#'   Test. Default is TRUE. Only used in the 2 by 2 by K case.
#' @param conf.level confidence level of the returned confidence 
#'   interval. Must be a single number between 0 and 1. Default is 0.95. Only used in the 2 by 2 by K case.
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'   \item{statistic}{The chi-square test statistic.} 
#'   \item{parameter}{The degrees of freedom of the approximate chi-squared distribution of the test statistic.}
#'   \item{p.value}{The p-value of the test.} 
#'   \item{estimate}{Common hazard ratio estimate. Only present in the 2 by 2 by K case.} 
#'   \item{null.value}{The null common hazard ratio, which is currently set to 1. Only present in the 2 by 2 by K case.} 
#'   \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided". Only present in the 2 by 2 by K case.} 
#'   \item{method}{A character string indicating the method employed.} 
#'   \item{data.name}{A character string giving the name of the data.} }
#' @seealso \code{\link{hazard.ratio.test}} for unstratified survival data and \code{\link{hazard.ratios}} 
#' for calculating estimated hazard ratios of stratified survival data.
#' @references Newman (2001), page 219 - 221, 226.
#' @examples 
#' ## Example 10.14
#' with(breast.survival, 
#'      k.hazard.ratio.test(time = time, status = status, 
#'                          exposure = receptor.level, strata = stage))
#'                          
#' ## LR Test 
#' with(breast.survival, 
#'      k.hazard.ratio.test(time = time, status = status,
#'                          exposure = receptor.level, strata = stage, 
#'                          Wald=FALSE))
#'                          
#' ## Example 10.18 - polychotomous exposure
#' with(breast.survival, 
#'      k.hazard.ratio.test(time = time, status = status,
#'                          exposure = stage, strata = receptor.level))
k.hazard.ratio.test <- function(time, status, exposure, strata, 
                                Wald=TRUE, conf.level=0.95){
  dname <- paste("\n  ", deparse(substitute(time)), 
                 "\n  ", deparse(substitute(status)))
  alternative <- "two.sided"
  
  # 2 exposure levels
  if(length(levels(exposure))==2){
    
    dk <- tapply(status, list(exposure, strata), sum)
    nk <- tapply(time, list(exposure, strata), sum)
    mk <- apply(dk,2,sum)
    n <- apply(nk,2,sum)
    
    # estimate hazard ratio
    f1 <- function(x){
      sum((x*mk*nk[1,])/(x*nk[1,] + nk[2,])) - sum(dk[1,])
    }
    est <- uniroot(f = f1, interval = c(0,1e5))$root
    names(est) <- "common hazard ratio"
    null <- 1
    names(null) <- names(est)
    
    hr2 <- mk/(est*nk[1,] + nk[2,])
    hr1 <- est*hr2
    # fitted counts:
    dh1 <- hr1*nk[1,]
    dh2 <- hr2*nk[2,]
    
    # confidence interval
    V <- sum((1/dh1 + 1/dh2)^(-1))
    alpha <- (1-conf.level)/2
    CINT <- exp(log(est) + c(-1,1)*qnorm(1 - alpha)/sqrt(V))
    attr(CINT, "conf.level") <- conf.level
    
    # tests of association
    # Null: log(HR) = 0
    # wald test of association
    if(Wald){
      V0 <- sum((nk[1,]*nk[2,]*mk)/(n^2))
      STATISTIC <- log(est)^2 * V0
      p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
      
    } else {
      # LRT of association
      # expected counts
      e1 <- nk[1,]*mk/n
      e2 <- nk[2,]*mk/n
      STATISTIC <- 2 * sum(dk[1,] * log((dk[1,]/e1)) + dk[2,] * log((dk[2,]/e2)))
      p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
    }
    names(STATISTIC) <- "X-squared"
    METHOD <- paste(if(Wald) "Wald" else "Likelihood Ratio", "Test of association")
    RVAL <- list(statistic = STATISTIC, parameter = c(df = 1), p.value = p.value, estimate = est, 
                 null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = METHOD, 
                 data.name = dname)   
  } else {
    dik <- tapply(status, list(exposure, strata), sum)
    nik <- tapply(time, list(exposure, strata), sum)
    mk <- apply(dik,2,sum)
    nk <- apply(nik,2,sum)
    # expected counts
    eik <- sapply(seq_along(levels(strata)),function(x)nik[,x]*mk[x]/nk[x])
    
    di. <- apply(dik,1,sum)
    ei. <- apply(eik,1,sum)
    
    STATISTIC <- sum(((di. - ei.)^2)/ei.)
    df <- nrow(dik) - 1
    names(df) <- "df"
    p.value <- pchisq(STATISTIC, df = df, lower.tail = FALSE)
    names(STATISTIC) <- "X-squared"
    METHOD <- ("X-squared Test of association - polychotomous exposure")
    RVAL <- list(statistic = STATISTIC, parameter = df, p.value = p.value, 
                 method = METHOD, 
                 data.name = dname)   
  }
  
  class(RVAL) <- "htest"
  return(RVAL)  
}


