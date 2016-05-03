#' Standardized Rate Ratio
#' 
#' Calculate the ratio of standardized death rates and its confidence interval.
#' 
#' @param count vector of age-specific cohort deaths
#' @param pop vector of age-specific cohort population or person-years
#' @param stdcount vector of age-specific standard population deaths
#' @param stdpop vector of age-specific standard population 
#' @param conf.level confidence level as a number between 0 and 1. Default is 0.95.
#'
#' @return a list with the following components:
#' \describe{
#' \item{Standardized rate ratio}{The ratio of the standardized death rates}
#' \item{CI}{The confidence interval for the standardized rate ratio}
#' \item{Crude rate ratio}{The ratio of the crude cohort death rate to the standard population death rate}
#' \item{directly standardized death rate}{The weighted average of the cohort age-specific death rates}
#' }
#'
#' @references Newman (2001), pages 251-255.
#' @seealso \code{\link[epitools]{ageadjust.direct}} in the \code{epitools} package.
#' @examples
#' with(schizophrenia, 
#'      std.rate.ratio(count = cohort.deaths, 
#'                     pop = cohort.py,
#'                     stdcount = alberta.deaths, 
#'                     stdpop = alberta.pop))
#' 
std.rate.ratio <- function(count, pop, stdcount, stdpop, conf.level=0.95){
  Ra <- sum((pop/sum(pop)) * 
              (count/pop))
  Rb <- sum((stdpop/sum(stdpop)) * 
              (stdcount/stdpop))
  
  # directly standardized death rate
  Ras <- sum((stdpop/sum(stdpop)) * 
               (count/pop))
  
  # crude rate ratio
  CRR <- Ra/Rb
  
  # standardized rate ratio
  SRR <- Ras/Rb
  
  # var of standardized death rate
  Nsk <- stdpop
  Ns <- sum(stdpop)
  Dak <- count
  Nak <- pop
  
  Dbk <- stdcount
  Nbk <- stdpop
  
  # sqrt(sum((Nsk/Ns)^2 * (Dak/(Nak^2))))
  varRas <- sum((Nsk/Ns)^2 * (Dak/(Nak^2)))
  varRbs <- sum((Nsk/Ns)^2 * (Dbk/(Nbk^2)))
  varlogSRR <- varRas/(Ras^2) + varRbs/(Rb^2)
  
  # 95 CI for SRR
  alpha <- (1 - conf.level)/2
  CINT <- exp(log(SRR) + c(-1,1) * qnorm(1 - alpha) * sqrt(varlogSRR))
  list("Standardized rate ratio" = SRR, "CI" = CINT,
       "Crude rate ratio" = CRR, 
       "directly standardized death rate" = Ras)
}


