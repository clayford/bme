#' Power calculation for Hazard Ratio Test
#' 
#' Compute the power of the Hazard Ratio Test, or determine parameters to obtain a target power.
#' 
#' Exactly one of the parameters n and power must be passed as NULL, and 
#' that parameter is determined from the others.
#'
#' @param n number of subjects needed for study
#' @param hr Hazard Ratio (ie, the "ratio worth detecting")
#' @param p1 proportion of cohort with history of exposure
#' @param pi2 probability a member of unexposed cohort will die during follow-up
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param alternative one- or two-sided test. Can be abbreviated.
#' @param tol numerical tolerance used in root finding, the default providing (at least) 
#' four significant digits. Root finding refers to \code{\link[stats]{uniroot}}, which is used to find sample size given power.
#'
#' @return Object of class "power.htest", a list of the arguments (including the computed one) 
#' augmented with method and note elements.
#' @export
#'
#' @references Newman (2001), pages 286 - 287
#' @seealso \code{\link{hazard.ratio.test}}
#' @examples
#' ## Example 14.4
#' power.hr.test(hr = 2, p1 = 0.251, pi2 = 0.174, power = 0.8)
#' power.hr.test(n = 414, hr = 2, p1 = 0.251, pi2 = 0.174)
power.hr.test <- function(n = NULL, hr, p1, pi2, sig.level = 0.05, power = NULL,
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
  p2 <- 1 - p1
  p.body <- quote({
    pnorm(sqrt((p1*p2*log(HR)^2) * (p1 * (1 - (1 - pi2)^HR) + p2*pi2) * n) - 
            qnorm(sig.level/tside, lower.tail = FALSE))
  })
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(n)) 
    n <- uniroot(function(n) eval(p.body) - power, c(1, 1e+07), 
                 tol = tol, extendInt = "upX")$root
  NOTE <- "n is the number of subjects needed for the study;\n      m is the number of deaths needed for the study"  
  METHOD <- "Hazard Ratio Test power calculation"
  structure(list(n = n, hr = hr, 
                 m = (qnorm(sig.level/tside, lower.tail = FALSE) + qnorm(power))^2 /(p1*p2*log(HR)^2),
                 sig.level = sig.level, 
                 power = power, alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}

