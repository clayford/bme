#' Mantel-Haenszel Test of association for matched-pairs case-control
#' 
#' Performs the Mantel-Haenszel Test of association for matched-pairs case-control data. 
#'
#' @param data a matrix with two rows and at least two columns with exposed cases on the 
#' first row and exposed controls on the first column (see data used in examples below)
#' @param conf.level confidence level of the returned confidence 
#'   interval. Must be a single number between 0 and 1. Default is 0.95.
#'
#' @return  A list with class \code{"htest"} containing the following 
#'   components: 
#'   \describe{ 
#'    \item{statistic}{The Mantel-Haenszel chi-square test statistic.} 
#'    \item{parameter}{The degrees of freedom of the test statistic.}
#'    \item{p.value}{The p-value of the test.} 
#'    \item{estimate}{The Mantel-Haenszel odds ratio estimate.} 
#'    \item{null.value}{The null hazard ratio, which is currently set to 1.} 
#'    \item{alternative}{A character string describing the alternative hypothesis. Currently only "two.sided".} 
#'    \item{method}{A character string indicating the method employed.} 
#'    \item{data.name}{A character string giving the name of the data.} 
#'   }
#' @references Newman (2001), pages 243 - 246.
#' @seealso \code{\link[stats]{mcnemar.test}}
#' @examples
#' ## Example 11.3 - (1:1) matched case-control data
#' estrogen
#' mantelhaen.mpcc.test(estrogen)
#' 
#' ## Example 11.4 - (1:M) matched case-control data
#' estrogen2
#' mantelhaen.mpcc.test(estrogen2)
mantelhaen.mpcc.test <- function(data, conf.level=0.95){
  dname <- deparse(substitute(data))
  alternative <- "two.sided"
  M <- ncol(data) - 1
  
  if(M==1){
    est <- data[1,2]/data[2,1]
    names(est) <- "odds ratio"
    null <- 1
    names(null) <- names(est)
    alpha <- (1-conf.level)/2
    v <- 1/data[1,2] + 1/data[2,1]
    CINT <- exp(log(est) + c(-1,1)*qnorm(1 - alpha)*sqrt(v))
    attr(CINT, "conf.level") <- conf.level
    # test of association
    STATISTIC <- ((data[1,2] - data[2,1])^2)/(data[1,2] + data[2,1])
    p.value <- pchisq(q = STATISTIC, df = 1, lower.tail = FALSE)
    
    names(STATISTIC) <- "X-squared"
    
  } else {
    # (1:M) matched case-control data
    R <- (1/(M+1))*sum(data[1,(1:M)]*(M + 1 - seq(M)))
    S <- (1/(M+1))*sum(data[2,(2:(M+1))]*seq(M))
    est <- R/S
    names(est) <- "odds ratio"
    null <- 1
    names(null) <- names(est)
    alpha <- (1-conf.level)/2
    T <- (1/(M+1)^2)*sum(data[1,(1:M)] * (M + 1 - seq(M)) * (M + 2 - seq(M)))
    U <- (1/(M+1)^2)*sum(data[2,(2:(M+1))]*seq(M) * (M - seq(M)))
    V <- (1/(M+1)^2)*sum(data[1,(1:M)] * (seq(M) - 1) * (M + 1 - seq(M)))
    W <- (1/(M+1)^2)*sum(data[2,(2:(M+1))] * seq(M) * (seq(M) + 1))
    # RBG variance estimate
    var_log_OR_mh <- T/(2*R^2) + (U + V)/(2*R*S) + W/(2*S^2)
    
    # 95% CI
    CINT <- exp(log(est) + c(-1,1)*qnorm(1 - alpha)*sqrt(var_log_OR_mh))
    
    # MH test of association
    num <- (sum(data[1,(1:M)]) - 
              sum((data[1,(1:M)] + data[2,(2:(M+1))]) * seq(M) / (M + 1)))^2
    
    den <- sum((data[1,(1:M)] + data[2,(2:(M+1))]) * seq(M) * (M + 1 - seq(M)) / (M + 1)^2)
    STATISTIC <- num/den
    p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
    names(STATISTIC) <- "X-squared"
  }
  METHOD <- paste("Mantel-Haenszel Test of association for matched-pairs case-control")
  RVAL <- list(statistic = STATISTIC, parameter = c(df = 1), p.value = p.value, 
               estimate = est, null.value = null,
               conf.int = CINT, alternative = alternative,
               method = METHOD, 
               data.name = dname) 
  class(RVAL) <- "htest"
  return(RVAL)
  
}

