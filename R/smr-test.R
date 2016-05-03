#' Standardized Mortality Ratio (SMR) Test
#' 
#' Performs a chi-sqaure test of the null that there is no mortality difference 
#' between the cohort and standard population.
#' 
#' When the expected number of deaths is less than 5, \code{\link{exact.rate.test}} 
#' is used to calculate the confidence interval and p-value. 
#'  
#'
#' @param count vector of age-specific cohort deaths
#' @param pop vector of age-specific cohort population or person-years
#' @param stdcount vector of age-specific standard population deaths
#' @param stdpop vector of age-specific standard population 
#' @param conf.level confidence level as a number between 0 and 1. Default is 0.95.
#'
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'    \item{p.value}{The p-value of the test.} 
#'    \item{estimate}{An estimate of the SMR ratio.} 
#'    \item{null.value}{The null SMR ratio. Currently set to 1.} 
#'    \item{conf.int}{A confidence interval for the SMR ratio.} 
#'    \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided".} 
#'    \item{method}{A character string indicating the method employed.} 
#'    \item{data.name}{Character strings giving the name of the data.} 
#'   }
#' @references Newman (2001), page 255-258.
#' @seealso \code{\link[epitools]{ageadjust.indirect}} in the \code{epitools} package.
#'
#' @examples
#' ## Example 12.3
#' 
#' with(schizophrenia, 
#'      smr.test(count = cohort.deaths,
#'               pop = cohort.py, 
#'               stdcount = alberta.deaths, 
#'               stdpop = alberta.pop))
#'                         
#' ## Example 12.3 continued,
#' ## Exact test (expected number of deaths < 5)
#' 
#' with(subset(schizophrenia,age.group=="10-19"), 
#'      smr.test(count = cohort.deaths,
#'               pop = cohort.py, 
#'               stdcount = alberta.deaths, 
#'               stdpop = alberta.pop))
#' 
smr.test <- function(count, pop, stdcount, stdpop, conf.level = 0.95){
  dname <- paste("\n  ", deparse(substitute(count)), 
                 "\n  ", deparse(substitute(pop)),
                 "\n  ", deparse(substitute(stdcount)),
                 "\n  ", deparse(substitute(stdpop)))
  alternative <- "two.sided"
  
  Da <- sum(count)
  Rsk <- stdcount/stdpop
  Nak <- pop
  Ea <- sum(Rsk*Nak)
  # standardized mortality ratio
  est <- Da/Ea
  names(est) <- "Standardized Mortality Ratio"
  null <- 1
  names(null) <- names(est)
  
  if(Ea < 5){
    CINT <- exact.rate.test(time = Ea, status = Da)$conf.int
    attr(CINT, "conf.level") <- conf.level
    p.value <- exact.rate.test(time = Ea, status = Da)$p.value
    METHOD <- paste("Exact test of no mortality difference between cohort and standard population")
    RVAL <- list(p.value = p.value, 
                 estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = METHOD, 
                 data.name = dname) 
    
  } else {
    varSMR <- est/Ea
    # 95% CI
    alpha <- (1-conf.level)/2
    CINT <- est + c(-1,1)*qnorm(1 - alpha)*sqrt(varSMR)
    attr(CINT, "conf.level") <- conf.level
    STATISTIC <- ((Da - Ea)^2)/Ea
    p.value <- pchisq(q = STATISTIC, df = 1, lower.tail = FALSE)
    names(STATISTIC) <- "X-squared"
    METHOD <- paste("Test of no mortality difference between cohort and standard population")
    RVAL <- list(statistic = STATISTIC, parameter = c(df = 1), p.value = p.value, 
                 estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = METHOD, 
                 data.name = dname) 
    
  }
  class(RVAL) <- "htest"
  return(RVAL)
}

