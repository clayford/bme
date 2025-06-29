#' Power calculations for unstratified closed cohort studies
#' 
#' Compute the power of risk difference, risk ratio and odds ratio methods for closed cohort studies, 
#' or determine parameters to obtain a target power.
#' 
#' Exactly one of the parameters n and power must be passed as NULL, and that parameter is determined 
#' from the others. 
#'
#' @param n number of observations
#' @param delta The true risk difference (RD), risk ratio (RR) or odds ratio (OR), depending on \code{type} selected. 
#' @param p2 probability that someone without a history of exposure will develop the disease
#' @param rho The ratio of unexposed to exposed subjects
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param type type of measurement: \code{RD} = "Risk Difference", \code{RR} = "Risk Ratio", and \code{OR} = "Odds Ratio"
#' @param alternative one- or two-sided test. Can be abbreviated.
#' @param tol numerical tolerance used in root finding, the default providing (at least) four significant digits. 
#' Root finding refers to \code{\link[stats]{uniroot}}, which is used to find sample size given power.
#' 
#'
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented 
#' with method and note elements. Sample size is returned as \code{r1} and \code{r2}. r1 is number of exposed 
#' subjects and r2 is number of unexposed subjects. r2 is equal to r1 * rho. 
#' @export
#'
#' @references Newman (2001), pages 283 - 285
#' @seealso \code{\link[stats]{power.prop.test}} which can perform Risk Difference power calculations assuming \code{rho = 1}
#' @examples
#' ## Example 14.2 (row 1 of Table 14.2)
#' power.closed.cohort(delta = 0.01, p2 = 0.05, rho = 1, type = "RD", power = 0.80)
#' 
#' ## Table 14.2
#' sapply(c(0.01, 0.05, 0.10, 0.20, 0.30), 
#'        function(x)power.closed.cohort(delta = x, p2 = 0.05, rho = 1, power = 0.8)$r1)
#'        
#' ## Table 14.3
#' sapply(lapply(c(1:5,10,20), 
#'        function(x)power.closed.cohort(delta = 0.05, rho = x, power = 0.8, p2 = 0.05)),
#' function(x)x[1:2])
#' 
#' ## Example 14.5 (row 2 of Table 14.4)
#' power.closed.cohort(delta = 3, p2 = 0.05, rho = 1, power = 0.8, type = "OR")
#' 
#' ## Table 14.4
#' sapply(c(2:5,10), function(x)power.closed.cohort(delta = x, p2 = 0.05, 
#'                                                  rho = 1, power = 0.8, type = "OR")$r1)
#' 
#' ## Example 14.7
#' power.closed.cohort(n = 100, delta = 2, p2 = 0.10, rho = 2, type = "OR")
#' 
#' ## Example 14.8
#' power.closed.cohort(rho = 6.73, p2 = 0.048, delta = 3, power = 0.8, type = "OR")
power.closed.cohort <- function(n = NULL, delta, p2, rho = 1,
                                sig.level = 0.05, power = NULL, 
                                type = c("RD","RR","OR"),
                                alternative = c("two.sided", "one.sided"),  
                                tol = .Machine$double.eps^0.25)
{
  if (sum(sapply(list(n, power), is.null)) != 
      1) 
    stop("exactly one of 'n' and 'power' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || 
      any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  type <- match.arg(type)
  if(type == "RD") p1 <- p2 + delta
  else if(type == "RR") p1 <- p2 * delta
  else p1 <- delta*p2 / (delta*p2 + (1 - p2))
  p0 <- (p1 + (p2 * rho)) / (1 + rho)
  p.body <- quote({
    stats::pnorm((abs(p1 - p2) * sqrt(n) - (stats::qnorm(sig.level/tside, lower.tail = FALSE)) * 
             sqrt(p0 * (1 - p0) * ((1 + rho)/rho))) / 
            sqrt(p1 * (1 - p1) + ((p2 * (1 - p2))/rho)))
  })
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(n)) 
    n <- stats::uniroot(function(n) eval(p.body) - power, c(1, 1e+07), 
                 tol = tol, extendInt = "upX")$root
  NOTE <- "r1 is number of exposed subjects; r2 is number of unexposed subjects"  
  METHOD <- paste(switch(type, RD = "Risk Difference", RR = "Risk Ratio", OR = "Odds Ratio"), 
                  "power calculation")
  structure(list(r1 = round(n), r2 = round(n) * rho, type = type, delta = delta,
                 p2 = p2, sig.level = sig.level, 
                 power = power, alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}

