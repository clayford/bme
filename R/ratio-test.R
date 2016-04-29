#' Hazard Ratio Test of association for unstratified survival data
#' 
#' Performs a test of the null that the hazard ratio between two groups is 1.
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param exposure a factor vector with two levels indicating exposure. The
#'   first level is assumed to be the exposed condition.
#' @param Wald a logical indicating whether to use the Wald or Likelihood Ratio 
#'   Test. Default is TRUE.
#' @param conf.level confidence level of the returned confidence 
#'   interval. Must be a single number between 0 and 1. Default is 0.95.
#' @param rev a logical indicating to whether to reverse the order of the factor
#'   levels when estimating the hazard ratio.
#' @param exact a logical indicating whether the asymptotic (unconditional) test
#'   or the exact test should be computed. Default is FALSE. The \code{"Wald"}
#'   argument is ignored when set to TRUE.
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'    \item{statistic}{The chi-square test statistic. Only returned if \code{"exact = FALSE"}.} 
#'    \item{p.value}{The p-value of the test.} 
#'    \item{estimate}{Hazard ratio estimate.} 
#'    \item{null.value}{The null hazard ratio, which is currently set to 1.} 
#'    \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided".} 
#'    \item{method}{A character string indicating the method employed.} 
#'    \item{data.name}{A character string giving the name of the data.} 
#'   }
#' @references Newman (2001), page 206 - 213.
#' @examples 
#' ## Example 10.8
#' with(breast.survival, hazard.ratio.test(time, status, receptor.level))
#' with(breast.survival, hazard.ratio.test(time, status, receptor.level, Wald = FALSE))
#' ## Example 10.9
#' ## Receptor Level-Breast Cancer: Stage III
#' dat <- subset(breast.survival, stage=="III")
#' hazard.ratio.test(time = dat$time, status = dat$status, exposure = dat$receptor.level)
#' hazard.ratio.test(time = dat$time, status = dat$status, exposure = dat$receptor.level, 
#'                   Wald = FALSE)
#' ## Example 10.10
#' hazard.ratio.test(time = dat$time, status = dat$status, exposure = dat$receptor.level, 
#'                   exact = TRUE)

hazard.ratio.test <- function(time, status, exposure, Wald=TRUE, conf.level=0.95, 
                              exact=FALSE){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  n <- tapply(time, exposure, sum)
  d <- tapply(status, exposure, sum)
  if(d[1] == 0 || d[2] == 0) est <- (d[1] + 0.5) * n2 / (d[2]  + 0.5) * n[1]
  else est <- (d[1] * n[2]) / (d[2] * n[1])
  names(est) <- "hazard ratio"
  null <- 1
  names(null) <- names(est)
  alpha <- (1-conf.level)/2
  
  if(!exact){
    # hazard ratio
    # HR conf interval
    CINT <- exp(log(est) + c(-1,1)*qnorm(1 - alpha)*sqrt((1/d[1]) + (1/d[2])))
    attr(CINT, "conf.level") <- conf.level
    
    # Wald and LRT tests of association
    m <- est*(d[2]/n[2])*n[1] + (d[2]/n[2])*n[2]
    e1 <- (n[1]*m)/(n[1] + n[2])
    e2 <- (n[2]*m)/(n[1] + n[2])
    
    # Wald test
    if(Wald){
      STATISTIC <- (log(est)^2 * n[1] * n[2] * m)/((n[1] + n[2])^2)
      p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)  
    } else {
      # LRT
      STATISTIC <- 2*(d[1] * log(d[1]/e1) + d[2] * log(d[2]/e2))
      p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
    }
    names(STATISTIC) <- "X-squared"
    METHOD <- paste(if(Wald) "Wald" else "Likelihood Ratio", "Test of association")
    RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = METHOD, 
                 data.name = dname)  
  } else {
    m <- sum(d)
    f1 <- function(x)1 - pbinom(d[1]-1, size = m, prob = x) - alpha
    f2 <- function(x)pbinom(d[1], size = m, prob = x) - alpha
    f1.out <- uniroot(f1, interval = c(0,1))
    f2.out <- uniroot(f2, interval = c(0,1))
    CINTpi <- c(f1.out$root, f2.out$root)
    CINT <- CINTpi*n[2]/((1 - CINTpi)*n[1]) # HR CI
    attr(CINT, "conf.level") <- conf.level
    
    # test
    # exact test of association
    pi_0 <- (null * n[1]) / (null * n[1] + n[2])
    p.value <- min( min(pbinom(q = d[1], size = m, prob = pi_0), 
                        1 - pbinom(q = d[1]-1, size = m, prob = pi_0)) * 2, 1)
    RVAL <- list(p.value = p.value, estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = "Exact Hazard Ratio Test of Association", 
                 data.name = dname)
  }
  
  class(RVAL) <- "htest"
  return(RVAL)
}
