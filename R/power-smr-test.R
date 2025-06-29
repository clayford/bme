#' Power calculation for Standardized Mortality Ratio (SMR) Test
#' 
#' Compute the power of the Standardized Mortality Ratio (SMR) Test, or determine parameters to 
#' obtain a target power.
#' 
#' Exactly one of the parameters n and power must be passed as NULL, and that parameter is determined 
#' from the others.
#'
#' @param n n is amount of person-time needed for study
#' @param smr Stanardized Mortality Ratio (ie, the "ratio worth detecting")
#' @param r death rate in the standard population
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param alternative one- or two-sided test. Can be abbreviated.
#' @param tol numerical tolerance used in root finding, the default providing (at least) four significant digits. 
#' Root finding refers to \code{\link[stats]{uniroot}}, which is used to find sample size given power.
#'
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented 
#' with method and note elements.
#' @export
#'
#' @references Newman (2001), pages 285 - 286. 
#' @seealso \code{\link{smr.test}}
#' @examples
#' ## Example 14.3
#' power.smr.test(smr = 1.5, r = 7283/957247, power = 0.8)
#' 
power.smr.test <- function(n = NULL, smr, r, sig.level = 0.05, power = NULL,
                           alternative = c("two.sided", "one.sided"),  
                           tol = .Machine$double.eps^0.25){
  if (sum(sapply(list(n, power), is.null)) != 
      1) 
    stop("exactly one of 'n' and 'power' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || 
      any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  p.body <- quote({
    stats::pnorm(sqrt(r * n * 4 * (sqrt(smr) - 1)^2) - stats::qnorm(sig.level/tside, lower.tail = FALSE))
  })
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(n)) 
    n <- stats::uniroot(function(n) eval(p.body) - power, c(1, 1e+07), 
                 tol = tol, extendInt = "upX")$root
  NOTE <- "n is amount of person-time needed for the study;\n      Ea is the expected number of deaths needed for the study"  
  METHOD <- "Standardized Mortality Ratio (SMR) Test power calculation"
  structure(list(n = n, smr = smr, r = r,
                 Ea = (stats::qnorm(sig.level/tside, lower.tail = FALSE) + stats::qnorm(power))^2 /(4 * (sqrt(smr) - 1)^2),
                 sig.level = sig.level, 
                 power = power, alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}
