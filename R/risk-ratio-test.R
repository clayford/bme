#' Risk Ratio Test for single 2 x 2 table
#' 
#' Performs a test of the null that the risk ratio between two groups in a 2 x 2 table is 1.
#' 
#' @param data a 2 x 2 matrix.
#' @param Wald a logical indicating whether to use the Wald or Likelihood Ratio 
#'   Test. Default is TRUE. 
#' @param conf.level confidence level of the returned confidence 
#'   interval. Must be a single number between 0 and 1. Default is 0.95.
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'   \item{statistic}{The chi-square test statistic.} 
#'   \item{parameter}{The degrees of freedom of the approximate chi-squared distribution of the test statistic.}
#'   \item{p.value}{The p-value of the test.} 
#'   \item{estimate}{The unconditional maximum likelihood estimate of the risk ratio.} 
#'   \item{null.value}{The null risk ratio, which is currently set to 1.} 
#'   \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided".} 
#'   \item{method}{A character string indicating the method employed.} 
#'   \item{data.name}{A character string giving the name of the data.} }
#' @seealso \code{\link[epitools]{riskratio}} in the epitools package.
#' @export
#' @references Newman (2001), pages 143-144.
#' @examples 
#' ## Example 6.1
#' risk.ratio.test(breast)
#' risk.ratio.test(breast, Wald = FALSE)
risk.ratio.test <- function(data, conf.level=0.95, Wald=TRUE){
  dname <- deparse(substitute(data))
  alternative <- "two.sided"
  
  a <- data[1,]
  b <- data[2,]
  r <- apply(data,2,sum)
  m <- apply(data,1,sum)
  # risk ratio
  if(a[1]==0 || a[2]==0){
    est <- ((a[1] + 0.5) *r[2])/((a[2] + 0.5)*r[1])
  } else {
    est <- (a[1]*r[2])/(a[2]*r[1])  
  }
  names(est) <- "risk ratio"
  null <- 1
  names(null) <- names(est)
  # variance
  v <- b[1]/(a[1]*r[1]) + b[2]/(a[2]*r[2])
  alpha <- (1-conf.level)/2
  CINT <- exp(log(est) + c(-1,1) * stats::qnorm(1 - alpha) * sqrt(v))
  attr(CINT, "conf.level") <- conf.level
  if(Wald){
    STATISTIC <- ((log(est)^2) * r[1] * r[2] * m[1]) / (sum(data) * m[2])
    p.value <- stats::pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  } else {
    # LRT of association
    STATISTIC <- 2*sum(data * log(data/epitools::expected(data)))
    p.value <- stats::pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  }
  names(STATISTIC) <- "X-squared"
  METHOD <- paste(if(Wald) "Wald" else "Likelihood Ratio", "Test of association for risk ratio")
  RVAL <- list(statistic = STATISTIC, parameter = c(df = 1), p.value = p.value, estimate = est, 
               null.value = null,
               conf.int = CINT, alternative = alternative,
               method = METHOD, 
               data.name = dname) 
  class(RVAL) <- "htest"
  return(RVAL)
}