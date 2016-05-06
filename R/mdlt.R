#' Multiple Decrement Life Table (MDLT) for a given year
#' 
#' Construct a multiple decrement life table for a given year. A multiple decrement life table 
#' (MDLT) describes the mortality experience of the group of individuals in the 
#' OLT cohort who are "due to die" of a particular cause of death. The MDLT is
#' an example of a competing risks model.
#' 
#' @param group The age group vector as integer. For example, 0, 1, 5, ... where 0 indicates birth to one year, 
#' 1 indicates age 1 to age 5, and so forth. See data used with example below.
#' @param all.deaths The number of deaths from all causes for each age group as a numeric vector.
#' @param cause.deaths The number of deaths from a particular cause for each age group as a numeric vector.
#' @param pop The population for each age group as a numeric vector.
#' @param radix The number of individuals in the OLT birth cohort. Typically set to a large number like 100,000.
#'   
#' @return a \code{data.frame} with \code{length(group)} rows and 7 columns:
#' \describe{
#' \item{group}{Age group}
#' \item{qjk}{crude conditional probability of dying in age group}
#' \item{lkxj}{number of survivors to age x}
#' \item{djk}{number of deaths in age group}
#' \item{Ljk}{number of person years in age group}
#' \item{Tkxj}{number of person years after age x}
#' \item{ekxj}{life expectancy at age x}
#' 
#' }
#' @export
#' @seealso \code{\link{olt}} and \code{\link{cdlt}}
#' @references Newman (2001), pages 270-273.
#' @examples
#' ## Example 13.3
#' ## create multiple decrement life table for neoplasms
#' with(males, mdlt(group = age.group, all.deaths = deaths.all, 
#'                  cause.deaths = neoplasm.deaths, pop = pop))
#' 
#' mdlt.out <- with(males,
#'                  mdlt(group = age.group, all.deaths = deaths.all,
#'                  cause.deaths = neoplasm.deaths, pop = pop))
#'                  
#' ## 27.17% of males born in Canada in 1991 predicted to die of neoplasm
#' mdlt.out[1,"Ljk"] / 1e5 * 100
#' 
#' ## For an individual due to die of neoplasm, the life expectancy 
#' ## at birth is 73.27 years
#' mdlt.out[1,"ekxj"]
mdlt <- function(group, all.deaths, cause.deaths, pop, radix=1e5){
  xj <- group
  
  # step 1
  Rj <- all.deaths/pop
  nj <- diff(group)
  qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
  Djk <- cause.deaths 
  Dj <- all.deaths
  qjk <- (Djk/Dj) * qj
  
  # step 2
  pj <- 1 - qj
  l0 <- radix
  lxj <- cumprod(c(l0,pj[-length(pj)]))
  djk <- qjk * lxj
  
  # step 3
  lkxj <- rev(cumsum(rev(djk)))
  
  # step 4
  dj <- qj * lxj
  Lj <- dj/Rj
  nj <- diff(xj)
  Ljk <- c(
    ((lkxj[-length(lkxj)] + lkxj[-1]) * nj)/2,
    (Djk[length(Djk)]/Dj[length(Dj)]) * Lj[length(Lj)]
  )
  
  # step 5
  Tkxj <- rev(cumsum(rev(Ljk)))
  
  # step 6
  ekxj <- Tkxj/lkxj
  
  data.frame(group = xj, qjk = round(qjk,5), lkxj = round(lkxj), 
             djk = round(djk), Ljk = round(Ljk),
             Tkxj = round(Tkxj), ekxj = round(ekxj,2))
}
