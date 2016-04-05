# Biostats in epidemiology
# ch 5 work


# 5.1 ---------------------------------------------------------------------

# asymptotic unconditional methods

# table 5.3
t53 <- array(data = c(2,10,5,50,
                      9,13,17,57,
                      12,2,9,6),
             dim = c(2,2,3), dimnames = list(survival=c("dead","alive"),
                                             ReceptorLevel = c("low","high"),
                                             stage=c("I","II","III")))
library(epitools)
# table 5.4
apply(t53, 3, oddsratio, method="wald")
apply(t53, 3, oddsratio, method="wald", verbose=TRUE)

# table 4.5(a)
t45a <- apply(t53,c(1,2),sum)
# reported on p. 126
# overall odds ratio for the cohort
# "crude" odds ratio estimate
oddsratio(t(t45a), method = "wald", rev="both")


# p. 122: summary estimate of OR
# unconditional maximum likelihood estimate
m1j <- c(7,26,21)
r1j <- c(12,22,14)
r2j <- c(55,74,15)
a1j <- c(2,9,12)

OR <- function(oru, m1j, r1j, r2j, a1j){
  (-1/(2*(oru - 1))) * sum((-1*((m1j + r1j)*oru - m1j + r2j)) + 
                           sqrt((-1*((m1j + r1j)*oru - m1j + r2j))^2 - 
                                  4*(oru - 1)*(oru*m1j*r1j))) - sum(a1j)
}

ORu <- uniroot(OR, interval = c(2,3), m1j = m1j, r1j = r1j, r2j = r2j, a1j = a1j)
ORu$root
# 2.51

# a-hat_1j
ahat_j <- function(oru, m1j, r1j, r2j){
(-1/(2*(oru - 1))) * ((-1*((m1j + r1j)*oru - m1j + r2j)) + 
                           sqrt((-1*((m1j + r1j)*oru - m1j + r2j))^2 - 
                                  4*(oru - 1)*(oru*m1j*r1j)))
}
ahat <- ahat_j(oru = ORu$root, m1j = m1j, r1j = r1j, r2j = r2j)
# p-hats at bottom of page 127
(m1j - ahat)/r2j


# how can we generalize this function?
oru <- -4
sqrt((-1*((m1j + r1j)*oru - m1j + r2j))^2 - 4*(oru - 1)*(oru*m1j*r1j))

# try a bunch of ORs
or <- seq(-100, 100)
tout <- sapply(or,function(x)OR(oru = x, m1j = m1j, r1j = r1j, r2j = r2j, a1j = a1j)) 
# drop NaN
tout2 <- tout[!is.nan(tout)]
or2 <- or[!is.nan(tout)]
max(or2[tout2 < 0])
min(or2[tout2 > 0])

# how to get these out of an array
m1j <- c(7,26,21)
r1j <- c(12,22,14)
r2j <- c(55,74,15)
a1j <- c(2,9,12)

m1j <- apply(t53,c(1,3),sum)[1,]
r1j <- apply(t53,c(2,3),sum)[1,]
r2j <- apply(t53,c(2,3),sum)[2,]
a1j <- t53[1,1,]

# tests of association
# expected counts
ec <- apply(t53,3,epitools::expected)
# Wald 
V0 <- sum(apply((1/ec),2,sum)^(-1))
wald <- log(2.51)^2 * V0
pchisq(wald, df = 1, lower.tail = FALSE)
# LRT
ac <- matrix(t53,ncol = 3)
lrt <- 2 * sum(ac * log(ac/ec))
pchisq(lrt, df = 1, lower.tail = FALSE)


# fitted values
ahat1 <- ahat_j(oru = ORu$root, m1j = m1j, r1j = r1j, r2j = r2j)
ahat2 <- m1j - ahat1
bhat1 <- r1j - ahat1
bhat2 <- r2j - ahat2
fc <- rbind(ahat1, bhat1, ahat2, bhat2)

# stratum specific odds ratios (does not match book's)
spor <- sapply(suppressWarnings(apply(t53,3,epitools::oddsratio)), function(x)x$measure[2,1])

# tests of homogeneity
# LRT
2 * sum(ac * log(ac/fc))


# function for unconditional maximum likelihood
or.u <- function(x){
  m1j <- apply(x,c(1,3),sum)[1,]
  r1j <- apply(x,c(2,3),sum)[1,]
  r2j <- apply(x,c(2,3),sum)[2,]
  a1j <- x[1,1,]
  OR <- function(oru, m1j, r1j, r2j, a1j){
    (-1/(2*(oru - 1))) * sum((-1*((m1j + r1j)*oru - m1j + r2j)) + 
                               sqrt((-1*((m1j + r1j)*oru - m1j + r2j))^2 - 
                                      4*(oru - 1)*(oru*m1j*r1j))) - sum(a1j)
  }
  or <- seq(-100, 100)
  suppressWarnings(tout <- sapply(or,function(x)OR(oru = x, m1j = m1j, r1j = r1j, r2j = r2j, a1j = a1j)))
  # drop NaNs
  tout2 <- tout[!is.nan(tout)]
  or2 <- or[!is.nan(tout)]
  lower <- max(or2[tout2 < 0])
  upper <- min(or2[tout2 > 0])
  
  ORu <- uniroot(OR, interval = c(lower,upper), m1j = m1j, r1j = r1j, r2j = r2j, a1j = a1j)
  ORu$root
  
  ahat_j <- function(oru, m1j, r1j, r2j){
    (-1/(2*(oru - 1))) * ((-1*((m1j + r1j)*oru - m1j + r2j)) + 
                            sqrt((-1*((m1j + r1j)*oru - m1j + r2j))^2 - 
                                   4*(oru - 1)*(oru*m1j*r1j)))
  }
  ahat1 <- ahat_j(oru = ORu$root, m1j = m1j, r1j = r1j, r2j = r2j)
  ahat2 <- m1j - ahat1
  bhat1 <- r1j - ahat1
  bhat2 <- r2j - ahat2
  fc <- rbind(ahat1, bhat1, ahat2, bhat2)
  
  # variance
  Vu <- 1/sum((1/ahat1 + 1/ahat2 + 1/bhat1 + 1/bhat2)^(-1))
  # CI
  logCI <- log(ORu$root) + c(-1,1)*c(qnorm(0.975)*sqrt(Vu))
  # Tests of association
  # expected counts
  ec <- apply(x,3,epitools::expected)
  # Wald 
  V0 <- sum(apply((1/ec),2,sum)^(-1))
  wald <- log(ORu$root)^2 * V0
  # LRT
  ac <- matrix(x,ncol = 3)
  lrt <- 2 * sum(ac * log(ac/ec))
  # Test of homogeneity
  lrth <- 2 * sum(ac * log(ac/fc))

  list(OR = ORu$root, fitted = list(ahat1, ahat2, bhat1, bhat2), variance = Vu,
       CI = exp(logCI), 
       Wald.statistic = wald, Wald.p.value = pchisq(wald, df = 1, lower.tail = FALSE),
       LRT.association.statistic = lrt, LRTa.p.value = pchisq(lrt, df = 1, lower.tail = FALSE),
       LRT.homogeneity.statistic = lrth, LRTh.p.value = pchisq(lrth, df = dim(x)[3]-1, lower.tail = FALSE))
}



or.u(t53)

# 5.2 ---------------------------------------------------------------------

# asymptotic conditional methods for J (2 X 2) Tables

# 5.3 ---------------------------------------------------------------------

mtout <- mantelhaen.test(t53, correct = FALSE)
str(mtout)

# MH OR manually
t53
Rj <- apply(t53,3,function(x)prod(diag(x))/sum(x))
Sj <- apply(t53,3,function(x)prod(x[2,1],x[1,2])/sum(x))
sum(Rj)/sum(Sj)

library(epiR)
rval <- epi.2by2(dat = t53, method = "cohort.count",
                 conf.level = 0.95, units = 100, homogeneity = "breslow.day",
                 outcome = "as.columns")
print(rval)


# 5.7 ---------------------------------------------------------------------

# example 5.5
# table 5.10
t510 <- array(data = c(2,10,9,13,12,2,
                      5,50,17,57,9,6),
             dim = c(2,3,2), dimnames = list(survival=c("dead","alive"),
                                             stage=c("I","II","III"),
                                             ReceptorLevel = c("low","high")))
# test of association
mantelhaen.test(t510)

# Recall: The Mantel-Haenszel adjusted measures of association are valid when
# the measures of association across the different strata are similar
# (homogenous), that is when the test of homogeneity of the odds (risk) ratios
# is not significant.

# table 5.11
# need to swap alive/dead to get MH estimates in table 5.11
mantelhaen.test(t510[c(2,1),c(1,2),], correct = FALSE)
mantelhaen.test(t510[c(2,1),c(1,3),], correct = FALSE)

# Table 5.12
expected(t510[,,1])
expected(t510[,,2])


tout <- lbl_test(as.table(t510), scores = list("stage" = c(1,2,3)))
statistic(tout)^2


apply(t510, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))

# another way
library(coin)
cmh_test(object = as.table(t510))

mh.est <- function(x){
  Rj <- apply(x,3,function(x)prod(diag(x))/sum(x))
  Sj <- apply(x,3,function(x)prod(x[2,1],x[1,2])/sum(x))
  sum(Rj)/sum(Sj)
}
mh.est(t510[,c(1,2),])


epi.2by2(t510[,c(1,2),])




# set up table for epiR 
t510df <- as.data.frame(as.table(t510))
t510b <- xtabs(Freq ~ stage + survival + ReceptorLevel, data=t510df)

mantelhaen.test(t510b[c(1,2),,])

apply(t510b[c(1,2),,],3,epi.2by2)
epi.2by2(dat = t510b[c(1,2),,], method = "cohort.count", 
         conf.level = 0.95, units = 100,  homogeneity = "breslow.day", 
         outcome = "as.columns")

# test for linear trend
# p. 141 result
eij <- apply(t510, 3, epitools::expected)
aij <- matrix(t510[1,,],ncol=2)
si <- c(1,2,3)
Uj <- sum(si %*% (aij - eij[c(1,3,5),]))


v1 <- ( sum(t510[2,,1]) / (sum(t510[,,1]) - 1)) * (sum(si^2 * eij[c(1,3,5),1]) - sum(si * eij[c(1,3,5),1])^2 / sum(eij[c(1,3,5),1]))
v2 <- ( sum(t510[2,,2]) / (sum(t510[,,2]) - 1)) * (sum(si^2 * eij[c(1,3,5),2]) - sum(si * eij[c(1,3,5),2])^2 / sum(eij[c(1,3,5),2]))

Vj <- sum(v1,v2)

X_t <- (Uj^2)/Vj # statistic
pchisq(X_t, df = 1, lower.tail = FALSE) # p-value

# function: test for linear trend
# disease on rows (yes in first row); exposure on columns; stratum is 3rd layer
tlt <- function(x, s = seq(dim(x)[2])){
  eij <- apply(x, 3, epitools::expected)
  aij <- matrix(x[1,,],ncol=dim(x)[3])
  si <- s
  i <- seq(dim(x)[2]) + (seq(dim(x)[2]) - 1)
  Uj <- sum(si %*% (aij - eij[i,]))
  
  
  v1 <- ( sum(x[2,,1]) / (sum(x[,,1]) - 1)) * (sum(si^2 * eij[i,1]) - sum(si * eij[i,1])^2 / sum(eij[i,1]))
  v2 <- ( sum(x[2,,2]) / (sum(x[,,2]) - 1)) * (sum(si^2 * eij[i,2]) - sum(si * eij[i,2])^2 / sum(eij[i,2]))
  
  Vj <- sum(v1,v2)
  X_t <- (Uj^2)/Vj
  list(statistic = X_t, p.value = pchisq(X_t, df = 1, lower.tail = FALSE))
}

tlt(t510)
