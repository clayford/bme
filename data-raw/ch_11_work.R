# Biostats in epidemiology
# ch 10 work

# Example 11.1
epitools::oddsratio(oral)

# Example 11.2
mantelhaen.test(oral.age, correct = FALSE)

# for Breslow-Day test of homogeneity
epiR::epi.2by2(dat = oral.age)


# 11.2 --------------------------------------------------------------------

# Example 11.3
# odds ratio

OR <- estrogen[1,2]/estrogen[2,1]
v <- 1/estrogen[1,2] + 1/estrogen[2,1]
CINT <- exp(log(OR) + c(-1,1)*qnorm(0.975)*sqrt(v))
# test of association
STATISTIC <- ((estrogen[1,2] - estrogen[2,1])^2)/(estrogen[1,2] + estrogen[2,1])
p.value <- pchisq(q = STATISTIC, df = 1, lower.tail = FALSE)

# put into a function
or.mpcc <- function(data, conf.level=0.95){
  dname <- deparse(substitute(data))
  alternative <- "two.sided"
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
  METHOD <- paste("Mantel-Haenszel Test of association for matched-pairs case-control")
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, null.value = null,
               conf.int = CINT, alternative = alternative,
               method = METHOD, 
               data.name = dname)  
class(RVAL) <- "htest"
return(RVAL)

}
or.mpcc(estrogen)


# mcnemar test
mcnemar.test(estrogen)


# 11.3 --------------------------------------------------------------------

estrogen2

# Example 11.4
# Mantel-Haenszel odds ratio estimate
M <- ncol(estrogen2) - 1
R <- (1/(M+1))*sum(estrogen2[1,(1:M)]*(M + 1 - seq(M)))
S <- (1/(M+1))*sum(estrogen2[2,(2:(M+1))]*seq(M))
est <- R/S

T <- (1/(M+1)^2)*sum(estrogen2[1,(1:M)] * (M + 1 - seq(M)) * (M + 2 - seq(M)))
U <- (1/(M+1)^2)*sum(estrogen2[2,(2:(M+1))]*seq(M) * (M - seq(M)))
V <- (1/(M+1)^2)*sum(estrogen2[1,(1:M)] * (seq(M) - 1) * (M + 1 - seq(M)))
W <- (1/(M+1)^2)*sum(estrogen2[2,(2:(M+1))] * seq(M) * (seq(M) + 1))
# RBG variance estimate
var_log_OR_mh <- T/(2*R^2) + (U + V)/(2*R*S) + W/(2*S^2)

# 95% CI
CINT <- exp(log(est) + c(-1,1)*qnorm(0.975)*sqrt(var_log_OR_mh))

# MH test of association
num <- (sum(estrogen2[1,(1:M)]) - 
  sum((estrogen2[1,(1:M)] + estrogen2[2,(2:(M+1))]) * seq(M) / (M + 1)))^2

den <- sum((estrogen2[1,(1:M)] + estrogen2[2,(2:(M+1))]) * seq(M) * (M + 1 - seq(M)) / (M + 1)^2)
STATISTIC <- num/den
p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)

# add to previous mh test function
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
# 1:1
mantelhaen.mpcc.test(estrogen)
# 1:M
mantelhaen.mpcc.test(estrogen2)
