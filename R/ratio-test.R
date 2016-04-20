#' Hazard Ratio Test of association for unstratified survival data
#' 
#' Performs a test of the null that the hazard ratio between two groups is 1.
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @param group a factor vector with at least two levels.
#' @param Wald a logical indicating whether to use the Wald or Likelihood Ratio 
#'   Test. Default is TRUE.
#' @param conf.level conf.level confidence level of the returned confidence 
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
#' with(breast.survival, hazard.ratio.test(time, status, receptor, rev = TRUE))
#' with(breast.survival, hazard.ratio.test(time, status, receptor, rev = TRUE, Wald = FALSE))
#' ## Example 10.9
#' ## Receptor Level-Breast Cancer: Stage III
#' dat <- breast.survival[breast.survival$stage=="III",]
#' hazard.ratio.test(time = dat$time, status = dat$status, group = dat$receptor, rev = TRUE)
#' hazard.ratio.test(time = dat$time, status = dat$status, group = dat$receptor, Wald = FALSE, rev = TRUE)
#' ## Example 10.10
#' hazard.ratio.test(time = dat$time, status = dat$status, group = dat$receptor, 
#'                   rev = TRUE, exact = TRUE)

hazard.ratio.test <- function(time, status, group, Wald=TRUE, conf.level=0.95, rev=FALSE,
                              exact=FALSE){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  g <- levels(group)
  # reverse order of groups?
  if(rev) g <- rev(g)
  if(!exact){
    n1 <- sum(time[group==g[1]])
    n2 <- sum(time[group==g[2]])
    d1 <- sum(status[group==g[1]])
    d2 <- sum(status[group==g[2]])
    
    # hazard ratio
    if(d1 == 0 || d2 == 0) est <- (d1 + 0.5) * n2 / (d2  + 0.5) * n1
    else est <- (d1 * n2) / (d2 * n1)
    names(est) <- "hazard ratio"
    null <- 1
    names(null) <- names(est)
    # HR conf interval
    alpha <- (1-conf.level)/2
    CINT <- exp(log(est) + c(-1,1)*qnorm(1 - alpha)*sqrt((1/d1) + (1/d2)))
    attr(CINT, "conf.level") <- conf.level
    
    # Wald and LRT tests of association
    m <- est*(d2/n2)*n1 + (d2/n2)*n2
    e1 <- (n1*m)/(n1 + n2)
    e2 <- (n2*m)/(n1 + n2)
    
    # Wald test
    if(Wald){
      STATISTIC <- (log(est)^2 * n1 * n2 * m)/((n1 + n2)^2)
      p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)  
    } else {
      # LRT
      STATISTIC <- 2*(d1 * log(d1/e1) + d2 * log(d2/e2))
      p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
    }
    names(STATISTIC) <- "X-squared"
    METHOD <- paste(if(Wald) "Wald" else "Likelihood Ratio", "Test of association")
    RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = METHOD, 
                 data.name = dname)  
  } else {
    d1 <- sum(status[group==g[1]])
    m <- sum(status)
    d2<- m - d1
    n1 <- sum(time[group==g[1]])
    n <- sum(time)
    n2 <- n - n1
    # estimate
    if(d1 == 0 || d2 == 0) est <- (d1 + 0.5) * n2 / (d2  + 0.5) * n1
    else est <- (d1 * n2) / (d2 * n1)
    names(est) <- "hazard ratio"
    null <- 1
    names(null) <- names(est)
    # CI
    alpha <- (1-conf.level)/2
    f1 <- function(x)1 - pbinom(d1-1, size = m, prob = x) - alpha
    f2 <- function(x)pbinom(d1, size = m, prob = x) - alpha
    f1.out <- uniroot(f1, interval = c(0,1))
    f2.out <- uniroot(f2, interval = c(0,1))
    CINTpi <- c(f1.out$root, f2.out$root)
    CINT <- CINTpi*n2/((1 - CINTpi)*n1) # HR CI
    attr(CINT, "conf.level") <- conf.level
    
    # test
    # exact test of association
    pi_0 <- (null * n1) / (null * n1 + n2)
    p.value <- min( min(pbinom(q = d1, size = m, prob = pi_0), 1 - pbinom(q = d1-1, size = m, prob = pi_0)) * 2, 1)
    RVAL <- list(p.value = p.value, estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = "Exact Hazard Ratio Test of Association", 
                 data.name = dname)
  }
  
  class(RVAL) <- "htest"
  return(RVAL)
}
