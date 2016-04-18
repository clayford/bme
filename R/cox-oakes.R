#' Cox-Oakes Test of Exponentiality
#' 
#' Performs a Cox-Oakes test of exponentiality of the null that the shape 
#' parameter of a Weibull distribution equals 1.
#' 
#' Under the null the Weibull distribution simplifies to an exponential 
#' distribution, implying a constant hazard function. Large values of the test 
#' statistic provide evidence against the exponential assumption of a constant 
#' hazard function. According to Newman (2001): "The correct interpretation is 
#' as follows: Given that we have decided to fit the data using using a Weibull 
#' model, not rejecting the null means there is no reason not to choose the 
#' exponential model (which is a type of Weibull model)." (p. 197)
#' 
#' @param time a numeric vector of survival times.
#' @param status a numeric vector of censoring indicators, with 0 = censored and
#'   1 = dead.
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'    \item{statistic}{The Cox-Oakes test statistic.} 
#'    \item{p.value}{The p-value of the test.} 
#'    \item{estimate}{An estimate of the hazard rate.} 
#'    \item{method}{A character string indicating the method employed.} 
#'    \item{data.name}{A character string giving the name of the data.} 
#'   }
#' @references Newman (2001), page 197.
#' @examples 
#' ## Example 10.1
#' cox.oakes(time = breast.survival$time, status = breast.survival$status)
#' ## Provides moderate evidence that exponential assumption may not be satisfied.
#' 
#' ## A graphical assessment can also be performed by plotting the estimated 
#' ## exponential survival curve and the Kaplan-Meier curve and deciding subjectively
#' ## whether the latter appears to be exponential in appearance.
#' 
#' ## Fig 10.2(a), p. 198
#' require(flexsurv)
#' fit.exp <- flexsurvreg(formula = Surv(time, status) ~ 1, data = breast.survival, dist="exp")
#' plot(fit.exp, ci=FALSE)
#' 
cox.oakes <- function(time,status) {
  dname <- deparse(substitute(time))
  d <- sum(status) # number of deaths in cohort
  n <- sum(time) # total amount of time cohort was under observation
  lam0.hat <- d/n # deaths per person-month
  names(lam0.hat) <- "hazard rate"
  si <- log(lam0.hat * time)
  ei.hat <- lam0.hat * time
  
  num <- (d + sum(si * (status - ei.hat)))^2
  denom <- d + sum((si^2)*ei.hat) - (sum(si*ei.hat))^2 / d
  X.co <- num/denom
  p.value <- pchisq(X.co, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = c(statistic = X.co), p.value = p.value,
               estimate = lam0.hat,
               method = "Cox-Oakes Test of Exponentiality", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

