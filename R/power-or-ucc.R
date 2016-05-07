#' Power calculations for Odds Ratio test in a unmatched case control study
#' 
#' Compute the power of an odds ratio test or determine the sample size to obtain a target power.
#' 
#' Exactly one of the parameters \code{n} and \code{power} must be passed as NULL, and that
#' parameter is determined from the others. 
#' 
#' Currently, only the two-sided alternative is calculated. 
#'
#' @param n number of cases
#' @param power power of test (1 minus Type II error probability)
#' @param rho ratio of controls to cases
#' @param p2 probability that a control has a history of exposure
#' @param OR odds ratio
#' @param sig.level significance level (Type I error probability)
#' @param tol numerical tolerance used in root finding, the default providing (at least) four significant digits.
#'
#' @return Object of class "power.htest", a list of the arguments (including the computed one) 
#' augmented with method and note elements.
#' @export
#' 
#' @references Newman (2001) pages 287, 293.
#'
#' @examples
#' ## Example 14.5 (Table 14.4)
#' sapply(c(2,3,4,5,10), function(x) power.icc(power = 0.8, rho = 1, p2 = 0.05, OR = x)$n)
#' 
#' ## Example 14.7
#' power.icc(n = 100, rho = 2, OR = 2, p2 = 0.10)
#' 
#' ## Example 14.8
#' power.icc(power = 0.8, rho = 6.73, OR = 3, p2 = 0.048)
power.or.ucc <- function (n = NULL, power = NULL, rho = 1, p2, OR, sig.level = 0.05,  
                          tol = .Machine$double.eps^0.25) 
{
  if (sum(sapply(list(n, power), is.null)) != 
      1) 
    stop("exactly one of 'n', and 'power'must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  alternative <- "two.sided"
  p1 <- (OR * p2)/(OR*p2 + (1 - p2))
  p0 <- (p1 + p2*rho)/(1 + rho)
  p.body <- quote({
    p1 <- (OR * p2)/(OR*p2 + (1 - p2))
    p0 <- (p1 + p2*rho)/(1 + rho)
    pnorm(((sqrt(n) * abs(p1 - p2)) - qnorm(1 - (sig.level/2))*sqrt(p0*(1 - p0)*((1 + rho)/rho))) /
            sqrt(p1*(1 - p1) + ((p2*(1 - p2))/rho)))
  })
  
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(n)) 
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07), 
                 tol = .Machine$double.eps^0.25, extendInt = "upX")$root
  else stop("internal error", domain = NA)
  NOTE <- "n is number of cases; multiply by rho to obtain number of controls"
  METHOD <- "Unmatched case-control power calculation"
  structure(list(n = n, power = power, p2 = p2, OR = OR, rho = rho, sig.level = sig.level, 
                 alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}



