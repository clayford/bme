#' Cause-Deleted Life Table (CDLT) for a given year
#' 
#' Construct a cause-deleted life table for a given year. The cause-deleted life table (CDLT)
#' is constructed just like an ordinary life table (OLT), except the cause-deleted death rate
#' is used in place of the annual death rate.
#'
#' @param group The age group vector as integer. For example, 0, 1, 5, ... where 0 indicates birth to one year, 
#' 1 indicates age 1 to age 5, and so forth. See data used with example below.
#' @param all.deaths The number of deaths from all causes for each age group as a numeric vector.
#' @param cause.deaths The number of deaths from a particular cause for each age group as a numeric vector.
#' @param pop The population for each age group as a numeric vector.
#' @param radix The number of individuals in the OLT birth cohort. Typically set to a large number like 100,000.
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
#' @seealso \code{\link{olt}} and \code{\link{mdlt}}
#' @references Newman (2001), pages 274-276.
#' @examples
#' ## Example 13.4
#' 
#' ## save ordinary life table
#' olt.out <- with(males, 
#'                 olt(group = age.group, deaths = deaths.all, pop = pop))
#'                 
#' ## save cause-deleted life table for circulatory disease
#' circ.cdlt <- with(males, 
#'                   cdlt(group = age.group, all.deaths = deaths.all, 
#'                        cause.deaths = circ.deaths, pop = pop))
#'  
#' ## save multiple decrement life table for circulatory disease
#' circ.mdlt <- with(males,
#'                   mdlt(group = age.group, all.deaths = deaths.all, 
#'                        cause.deaths = circ.deaths, pop = pop))
#'                        
#' ## Table 13.8 calculations for circulatory disease
#' 
#' ## Circulatory disease accounts for 40.09% of deaths
#' circ.mdlt[1,"lkxj"]/olt.out[1,"lxj"] * 100
#' 
#' ## life expectancy (at birth) for those due to die of circulatory disease 
#' ## is 77.7 years
#' circ.mdlt[1,"ekxj"]
#' 
#' ## Eliminating circulatory diseases as a cause of death would increase 
#' ## overall life expectancy by 6.06 years
#' circ.cdlt[1,"exj"] - olt.out[1,"exj"]
#' 
#' ## As a result of eliminating circulatory disease, the probability that 
#' ## a member of the MDLT cohort will survive to age 65 increases by about 13%
#' (circ.cdlt[circ.cdlt$group==65,"lxj"] - olt.out[olt.out$group==65,"lxj"]) / 
#'   circ.mdlt[1,"lkxj"] * 100
#'  
#' ## Eliminating circulatory diseases as a cause of death would increase 
#' ## the life expectancy of those due to die of this cause by 15.12 years
#' (olt.out[1,"lxj"] * (circ.cdlt[1,"exj"] - olt.out[1,"exj"])) / circ.mdlt[1,"lkxj"]
cdlt <- function(group, all.deaths, cause.deaths, pop, radix=1e5){
  Djk <- all.deaths - cause.deaths
  # annual death rate in population
  Rj <- (Djk/pop)
  # length of age group
  nj <- diff(group)
  # conditional probability of dying 
  qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
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


