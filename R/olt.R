#' Ordinary Life Table (OLT) for a given year
#' 
#' Construct an ordinary life table for a given year. The ordinary life table (OLT) provides a method of analyzing
#' mortality for all causes of death combined.
#'
#' @param group The age group vector as integer. For example, 0, 1, 5, ... where 0 indicates birth to one year, 
#' 1 indicates age 1 to age 5, and so forth. See data used with example below.
#' @param deaths The number of deaths from all causes for each age group as a numeric vector.
#' @param pop The population for each age group as a numeric vector.
#' @param radix The number of individuals in the OLT birth cohort. Typically set to a large number like 100,000.
#' @param phi The multiplier for sensitivity analysis. For example, set to 0.9 to examine mortality 
#' implications for a declining death rate.
#' @param linear logical indicating which functional form to use to calculate number of survivors to age x (lxj).
#'  Default is TRUE, which means use linear form. Otherwise use exponential.
#'
#' @return a \code{data.frame} with \code{length(group)} rows and 8 columns:
#' \describe{
#' \item{group}{Age group}
#' \item{qj}{conditional probability of dying }
#' \item{pj}{conditional probability of surviving}
#' \item{lxj}{number of survivors to age x}
#' \item{dj}{number of deaths in age group}
#' \item{Lj}{number of person years in age group}
#' \item{Txj}{Number of person years after age x}
#' \item{exj}{Life expectancy at age x}
#' }
#' 
#' @seealso \code{\link{cdlt}} and \code{\link{mdlt}}
#' @references Newman (2001), pages 264-269.
#' @examples
#' ## Example 13.1
#' with(males, olt(group = age.group, deaths = deaths.all, pop = pop))
#' ## predicted life expectancy birth is 73.34
#' 
#' ## with sensitivity analysis for declinind death rate:
#' with(males, olt(group = age.group, deaths = deaths.all, pop = pop, phi = 0.9))
#' ## predicted life expectancy birth is 75.65
olt <- function(group, deaths, pop, radix=1e5, phi=1, linear = TRUE){
  # annual death rate in population
  Rj <- (deaths/pop)
  # multiply by phi for sensitivity analysis
  Rj <- Rj * phi
  # length of age group
  nj <- diff(group)
  # conditional probability of dying
  if(linear){
    qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
    # exponential form
  } else {
    qj <- 1 - exp(-1 * nj * Rj)
  }
  # conditional probability of surviving
  pj <- 1 - qj
  # the radix: number of individuals in the OLT birth cohort
  l0 <- radix
  # number of survivors to age x
  lxj <- cumprod(c(l0,pj[-length(pj)]))
  # number of deaths in age group
  dj <- qj * lxj
  # number of person years in age group
  Lj <- dj/Rj
  # Number of person years after age x
  Txj <- sapply(1:length(Lj), function(x)sum(Lj[x:length(Lj)]))
  # Life expectancy at age x
  exj <- Txj/lxj
  # collect results into data frame and display
  data.frame(group = group, qj = round(qj,5), pj = round(pj, 5), 
             lxj = round(lxj), dj = round(dj), Lj = round(Lj), 
             Txj = round(Txj), exj = round(exj,2))
}


