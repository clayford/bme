#' Power Calculation for Matched-Pairs Case-Control Odds Ratio Test
#' 
#' Compute the power of the Matched-Pairs Case-Control Odds Ratio Test, or determine parameters to obtain a target power.
#' 
#' Exactly one of the parameters n and power must be passed as NULL, and that parameter is determined from the others.
#' 
#' When pre-study information on concordance is limited, Fleiss and Levin (1988) recommend using a relatively 
#' large value of \code{v}, such as 2.5.
#' 
#' Sample size estimates for M >=1 can be overestimated when the probability is small that a sample control 
#' has a history of exposure. 
#' 
#'
#' @param n number of matched pairs needed for study
#' @param OR Odds Ratio (ie, the "ratio worth detecting")
#' @param p2 probability a control has a history of exposure
#' @param v measure of concordance (of exposure) in the case and control samples. 
#' \code{v} = 1 if and only if the matching variables are not associated with exposure in the case and control samples.
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param M number of controls per case (default = 1)
#' @param alternative one- or two-sided test. Can be abbreviated.
#' @param tol numerical tolerance used in root finding, the default providing (at least) four significant digits. 
#' Root finding refers to \code{\link[stats]{uniroot}}, which is used to find sample size given power.
#' 
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with 
#' method and note elements.
#' @export
#'
#' @references Newman (2001), pages 288 - 291
#' 
#' Dupont, WD. (1988) Power calculations for matched case-control studies. \emph{Biometrics} 44, 1157-1168.
#' 
#' Fleiss, JL and Levin B. (1988) Sample size determination in studies with matched pairs. \emph{Journal of Clinical Epidemiology} 41, 727-730.
#'  
#' 
#' @examples
#' ## Example 14.6
#' power.mpcc(power = 0.80, OR = 3, p2 = 0.05, v = 4.82)
power.mpcc <- function(n = NULL, OR, p2, v, sig.level = 0.05, power = NULL, M = 1,
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
    p1 <- OR*p2 /(OR*p2 + (1 - p2))
    p0 <- p2*(1 - p2)*(OR +1) / (OR*p2+ (1 - p2))
    num <- 2*(v - 1)*p1*(1 - p1)
    den <- sqrt(1 + 4*(v - 1)*p1*(1-p1)) - 1
    pnorm((sqrt(p0 * n * ((M+1)/2*M) * (den/num) * (OR - 1)^2) - 
             qnorm(sig.level/tside, lower.tail = FALSE)*(OR+1)) / (2 * sqrt(OR)))
  })
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(n)) 
    n <- uniroot(function(n) eval(p.body) - power, c(1, 1e+07), 
                 tol = tol, extendInt = "upX")$root
  NOTE <- "n is the number of matched pairs needed for the study;\n      r is the number of discordant pairs needed for the study"  
  METHOD <- "Matched Pairs Case Control power calculation"
  structure(list(n = n,  
                 r = (qnorm(sig.level/tside, lower.tail = FALSE)*(OR + 1) + 2*qnorm(power)*sqrt(OR))^2 /(OR - 1)^2,
                 v = v, OR = OR, M = M, sig.level = sig.level, 
                 power = power, alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}
# Example 14.6

