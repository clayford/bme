# Biostats in epidemiology
# ch 7 work

#Risk difference methods for closed cohort data
library(epitools)

t45a <- matrix(c(23,25,31,113),ncol=2, dimnames = list(survival = c("dead","alive"), 
                                                       receptor.level=c("low","high")))
# Example 7.1
# the following produces 95% CI and Wald test of association
prop.test(x = t(t45a), correct = FALSE)
abs(diff(prop.test(x = t(t45a), correct = FALSE)$estimate)) # RR

# LRT test of association
X_lr <- 2*sum(t45a*log(t45a/epitools::expected(t45a)))
pchisq(X_lr, 1, lower.tail = FALSE)


# 7.3 ---------------------------------------------------------------------

t53 <- array(data = c(2,10,5,50,
                      9,13,17,57,
                      12,2,9,6),
             dim = c(2,2,3), dimnames = list(survival=c("dead","alive"),
                                             ReceptorLevel = c("low","high"),
                                             stage=c("I","II","III")))
# Example 7.3
# Mantel-Haenszel estimate of the risk difference
r1j <- apply(t53[,1,],2,sum)
r2j <- apply(t53[,2,],2,sum)
a1j <- t53[1,1,]
a2j <- t53[1,2,]
b1j <- t53[2,1,]
b2j <- t53[2,2,]
rj <- apply(t53,3,sum)

R <- sum((a1j * r2j) / rj)
S <- sum((a2j * r1j) / rj)
T <- sum((r1j * r2j) / rj)
RD_mh <- (R - S) / Tj

# Variance of M-H estimate of RD

U <- sum((r1j^2 * a2j - r2j^2 * a1j + r1j * r2j * (r2j - r1j) / 2) / rj^2)
V <- sum((a1j * b2j + a2j * b1j) / (2 * rj))
varRD_mh <- ((RD_mh * U) + V) / T^2

# 95% CI
RD_mh + c(-1,1) * qnorm(0.975) * sqrt(varRD_mh)

# function for MH est of RD with 95% CI
# disease on rows (yes in first row); exposure on columns; stratum is 3rd layer 
mh.rd <- function(x){
  if(length(dim(x)) != 3 || dim(x)[-3] != c(2,2)) stop("This function is for 2 x 2 x K tables")
  r1j <- apply(x[,1,],2,sum)
  r2j <- apply(x[,2,],2,sum)
  a1j <- x[1,1,]
  a2j <- x[1,2,]
  b1j <- x[2,1,]
  b2j <- x[2,2,]
  rj <- apply(x,3,sum)
  
  R <- sum((a1j * r2j) / rj)
  S <- sum((a2j * r1j) / rj)
  T <- sum((r1j * r2j) / rj)
  RD_mh <- (R - S) / Tj
  
  # Variance of M-H estimate of RD
  
  U <- sum((r1j^2 * a2j - r2j^2 * a1j + r1j * r2j * (r2j - r1j) / 2) / rj^2)
  V <- sum((a1j * b2j + a2j * b1j) / (2 * rj))
  varRD_mh <- ((RD_mh * U) + V) / T^2
  
  list(MH.risk.difference = RD_mh,
       CI = RD_mh + c(-1,1) * qnorm(0.975) * sqrt(varRD_mh))
}
mh.rd(t53)