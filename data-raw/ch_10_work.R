# Biostats in epidemiology
# ch 10 work


# 10.1 --------------------------------------------------------------------

# Weibull hazard function
time <- seq(0,4,0.01)
weibull.hazard <- function(alpha, lambda, t){
  if(alpha <= 0 || lambda <= 0) stop("alpha and lambda must be positive.")
  alpha * lambda * (lambda * t)^(alpha - 1)
}

# wh1 <- weibull.hazard(alpha = 3.0, lambda = 0, t = time)
wh1 <- weibull.hazard(alpha = 3.0, lambda = 1, t = time)
wh2 <- weibull.hazard(alpha = 1.5, lambda = 1, t = time)
wh3 <- weibull.hazard(alpha = 1.0, lambda = 1, t = time)
wh4 <- weibull.hazard(alpha = 0.5, lambda = 1, t = time)

# Fig 10.1(a)
plot(time, wh1, ylim = c(0,4), ylab = "Hazard", xlab = "Time", type="l")
lines(time, wh2, lty = 2)
lines(time, wh3, lty = 3)
lines(time, wh4, lty = 4)
legend("topright", legend = paste("alpha =", c(3.0,1.5,1.0,0.5)), lty= 1:4)

# Weibull survival function
weibull.survival <- function(alpha, lambda, t){
  if(alpha <= 0 || lambda <= 0) stop("alpha and lambda must be positive.")
  exp(-(lambda * t)^alpha)
}

ws1 <- weibull.survival(alpha = 3.0, lambda = 1, t = time)
ws2 <- weibull.survival(alpha = 1.5, lambda = 1, t = time)
ws3 <- weibull.survival(alpha = 1.0, lambda = 1, t = time)
ws4 <- weibull.survival(alpha = 0.5, lambda = 1, t = time)
# Fig 10.1(b)
plot(time, ws1, ylim = c(0,1), ylab = "Survival Probability", xlab = "Time", type="l")
lines(time, ws2, lty = 2)
lines(time, ws3, lty = 3)
lines(time, ws4, lty = 4)
legend("topright", legend = paste("alpha =", c(3.0,1.5,1.0,0.5)), lty= 1:4)

# Example 10.1

#   survreg's scale  =    1/(rweibull shape)
#   survreg's intercept = log(rweibull scale)

dat <- read.csv("BreastCancer.csv")
library(flexsurv)
fitw <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat, dist="weibull")
# Fig 10.2(b)
plot(fitw, ci=FALSE)
fitws <- summary(fitw)
fitws[[1]]$est

fite <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat, dist="exp")
# Fig 10.2(a)
plot(fite, ci=FALSE)
fites <- summary(fite)
fites[[1]]$est



library(survival)
m1 <- survreg(Surv(time, status) ~ 1, data=dat, dist = "weibull")
sm1 <- summary(m1)
1/sm1$scale # alpha-hat (shape)
1/exp(sm1$coefficients[1]) # lambda-hat (rate)

# cox-oakes test of exponentiality
# Null: alpha = 1, Weibull shape parameter is 1, hazard function is constant
d <- sum(dat$status) # number of deaths in cohort
n <- sum(dat$time) # total amount of time cohort was under observation
lam0.hat <- d/n # deaths per person-month
si <- log(lam0.hat * dat$time)
ei.hat <- lam0.hat * dat$time

num <- (d + sum(si * (dat$status - ei.hat)))^2
denom <- d + sum((si^2)*ei.hat) - (sum(si*ei.hat))^2 / d
X.co <- num/denom
pchisq(X.co, df = 1, lower.tail = FALSE)
# provides moderate evidence that exponential assumption may not be satisfied

# Funciton: cox-oakes test
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

cox.oakes(time = breast.survival$time, status = breast.survival$status)
c.out <- cox.oakes(time = breast.survival$time, status = breast.survival$status)
str(c.out)

# example 10.2
dat2 <- read.csv("ovarian.csv")

# Weibull model
m2 <- survreg(Surv(time, status) ~ 1, data=dat2, dist = "weibull", subset = grade=="High")
sm2 <- summary(m2)
1/sm2$scale # alpha-hat (shape)
1/exp(sm2$coefficients[1]) # lambda-hat (rate)

# Exponential model
m3 <- survreg(Surv(time, status) ~ 1, data=dat2, dist = "exponential", subset = grade=="High")
sm3 <- summary(m3)
1/exp(sm3$coefficients[1]) # lambda-hat (rate); alpha-hat (shape) is 1

# Fig 10.3(a)
library(flexsurv)
fite <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat2, dist="exp")
plot(fite, ci=FALSE)

# Fig 10.3(b)
fitw <- flexsurvreg(formula = Surv(time, status) ~ 1, data = dat2, dist="weibull")
plot(fitw, ci=FALSE)

cox.oakes(time = dat2$time[dat2$grade=="High"], status = dat2$status[dat2$grade=="High"])

# exact methods for single sample
# Example 10.3

# with n assumed constant (where n is total person time in study), d is a
# poisson random variable with parameter nu = lambda * n

# Null: lambda_0 = 0.4
d <- 2
n <- 10

ppois(q = 2, lambda = n*0.4) * 2

# Table 10.1
d <- 0:12
`D_=_d` <- dpois(d, lambda = n * 0.4)
`D_<=_d` <- ppois(d, lambda = n * 0.4)
`D_>=_d` <- ppois(d - 1, lambda = n * 0.4, lower.tail = FALSE)
cbind(d, `D_=_d`, `D_<=_d`, `D_>=_d`)


1 - ppois(q = 1, lambda = 10*0.024)
ppois(q = 2, lambda = 10*0.723)

f1 <- function(x)1 - ppois(q = 1, lambda = 10*x) - 0.025
f2 <- function(x)ppois(q = 2, lambda = 10*x) - 0.025
uniroot(f1, interval = c(0,10))
uniroot(f2, interval = c(0,10))

# exact.rate.ci <- function(d,n,ci=0.95){
#   alpha <- (1-ci)/2
#   f1 <- function(x)1 - ppois(q = d-1, lambda = n*x) - alpha
#   f2 <- function(x)ppois(q = d, lambda = n*x) - alpha
#   f1.out <- uniroot(f1, interval = c(0,n))
#   f2.out <- uniroot(f2, interval = c(0,n))
#   list(lower = f1.out$root, upper = f2.out$root)
# }
# exact.rate.ci(d = 2, n = 10)
# exact.rate.ci(d = 2, n = 10, ci=0.90)


exact.rate.test <- function(time, status, null=0, conf.level = 0.95){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  # estimate
  d <- sum(status, na.rm = TRUE)
  n <- sum(time, na.rm = TRUE)
  est <- d/n
  names(est) <- "hazard rate"
  names(null) <- names(est)
  # CI
  alpha <- (1-conf.level)/2
  f1 <- function(x)1 - ppois(q = d-1, lambda = n*x) - alpha
  f2 <- function(x)ppois(q = d, lambda = n*x) - alpha
  f1.out <- uniroot(f1, interval = c(0,n))
  f2.out <- uniroot(f2, interval = c(0,n))
  CINT <- c(f1.out$root, f2.out$root)
  attr(CINT, "conf.level") <- conf.level
  
  # test
  p.value <- min( min(ppois(q = d, lambda = n*null), 1 - ppois(q = d-1, lambda = n*null)) * 2, 1)
  RVAL <- list(p.value = p.value, estimate = est, null.value = null,
       conf.int = CINT, alternative = alternative,
       method = "Exact Hazard Rate Test for a single sample", 
       data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

exact.rate.test(time = 10, status = 2, null = 0.4)
with(breast.survival, exact.rate.test(time = time, status = status))
with(breast.survival, exact.rate.test(time = time, status = status, null = 0.01))
hout <- with(breast.survival, exact.rate.test(time = time, status = status, null = 0.01))
str(hout)


# asympotic methods for single sample
rate.test <- function(time, status, null=1, conf.level = 0.95, explicit=TRUE){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  # estimate
  d <- sum(status, na.rm = TRUE)
  n <- sum(time, na.rm = TRUE)
  est <- d/n
  names(est) <- "hazard rate"
  names(null) <- names(est)
  # CI
  alpha <- (1-conf.level)/2
  if(explicit){
    CINT <- est + c(-1,1)*((qnorm(1 - alpha)*sqrt(d))/n)  
  } else {
    a <- n^2
    b <- -n*(2*d + qnorm(1 - alpha)^2)
    c <- d^2
    CINT <- (-b + c(-1,1)*sqrt(b^2 - 4*a*c))/(2*a)
  }
  attr(CINT, "conf.level") <- conf.level
  
  # test
  STATISTIC <- ((d - null*n)^2)/(null*n)
  names(STATISTIC) <- "X-squared"
  p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, 
               null.value = null,
               conf.int = CINT, alternative = alternative,
               method = "Asymptotic Hazard Rate Test for a single sample", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}
# Examples 10.5 and 10.6
rate.test(time = 10, status = 2, null = 0.4)
rate.test(time = 25, status = 5, null = 0.4)
rate.test(time = 50, status = 10, null = 0.4)

rate.test(time = 10, status = 2, null = 0.4, explicit = FALSE)
rate.test(time = 25, status = 5, null = 0.4, explicit = FALSE)
rate.test(time = 50, status = 10, null = 0.4, explicit = FALSE)

# Example 10.7
with(breast.survival, rate.test(time = time, status = status))
with(breast.survival, rate.test(time = time, status = status, explicit = FALSE))



# 10.2 --------------------------------------------------------------------

# Poisson methods for unstratified survival data

# Example 10.8
n1 <- sum(breast.survival$time[breast.survival$receptor=="Low"])
n2 <- sum(breast.survival$time[breast.survival$receptor=="High"])
d1 <- sum(breast.survival$status[breast.survival$receptor=="Low"])
d2 <- sum(breast.survival$status[breast.survival$receptor=="High"])

# hazard ratio
hr.hat <- (d1 * n2)/(d2 * n1)

# cox-oakes for low and high receptor levels
with(breast.survival[breast.survival$receptor=="Low",], 
     cox.oakes(time = time, status = status))
with(breast.survival[breast.survival$receptor=="High",], 
     cox.oakes(time = time, status = status))

# stratum specific HR estimates
with(breast.survival[breast.survival$receptor=="Low",], 
     rate.test(time = time, status = status, explicit = FALSE))
with(breast.survival[breast.survival$receptor=="High",], 
     rate.test(time = time, status = status, explicit = FALSE))

# HR conf interval
exp(log(hr.hat) + c(-1,1)*qnorm(0.975)*sqrt((1/d1) + (1/d2)))

# Wald and LRT tests of association
# we say there is no association between exposure and survival if lamba_1 = lambda_2

m <- hr.hat*(d2/n2)*n1 + (d2/n2)*n2
e1 <- (n1*m)/(n1 + n2)
e2 <- (n2*m)/(n1 + n2)

# Wald test
X_w <- (log(hr.hat)^2 * n1 * n2 * m)/((n1 + n2)^2)
pchisq(X_w, df = 1, lower.tail = FALSE)
# LRT
X_lr <- 2*(d1 * log(d1/e1) + d2 * log(d2/e2))
pchisq(X_lr, df = 1, lower.tail = FALSE)

# Example 10.9
dat <- breast.survival[breast.survival$stage=="III",]
n1 <- sum(dat$time[dat$receptor=="Low"])
n2 <- sum(dat$time[dat$receptor=="High"])
d1 <- sum(dat$status[dat$receptor=="Low"])
d2 <- sum(dat$status[dat$receptor=="High"])

# hazard ratio
hr.hat <- (d1 * n2)/(d2 * n1)
# HR conf interval
exp(log(hr.hat) + c(-1,1)*qnorm(0.975)*sqrt((1/d1) + (1/d2)))

m <- hr.hat*(d2/n2)*n1 + (d2/n2)*n2
e1 <- (n1*m)/(n1 + n2)
e2 <- (n2*m)/(n1 + n2)

# Wald test
X_w <- (log(hr.hat)^2 * n1 * n2 * m)/((n1 + n2)^2)
pchisq(X_w, df = 1, lower.tail = FALSE)
# LRT
X_lr <- 2*(d1 * log(d1/e1) + d2 * log(d2/e2))
pchisq(X_lr, df = 1, lower.tail = FALSE)
# moderate evidence for an association between receptor level and survival in the stage II cohort


hazard.ratio.test <- function(time, status, group, Wald=TRUE, conf.level=0.95, rev=FALSE){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  g <- levels(group)
  # reverse order of groups?
  if(rev) g <- rev(g)
  n1 <- sum(time[group==g[1]])
  n2 <- sum(time[group==g[2]])
  d1 <- sum(status[group==g[1]])
  d2 <- sum(status[group==g[2]])
  
  # hazard ratio
  if(d1 == 0 || d2 == 0) est <- (d1 + 0.5) * n2 / (d2  + 0.5) * n1
  else est <- d1 * n2 / d2 * n1
  names(est) <- "hazard ratio"
  
    # HR conf interval
  alpha <- (1-conf.level)/2
  CINT <- exp(log(hr.hat) + c(-1,1)*qnorm(1 - alpha)*sqrt((1/d1) + (1/d2)))
  attr(CINT, "conf.level") <- conf.level
  
  # Wald and LRT tests of association
  m <- hr.hat*(d2/n2)*n1 + (d2/n2)*n2
  e1 <- (n1*m)/(n1 + n2)
  e2 <- (n2*m)/(n1 + n2)
  
  # Wald test
  if(Wald){
    STATISTIC <- (log(hr.hat)^2 * n1 * n2 * m)/((n1 + n2)^2)
    p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)  
  } else {
    # LRT
    STATISTIC <- 2*(d1 * log(d1/e1) + d2 * log(d2/e2))
    p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  }
  names(STATISTIC) <- "X-squared"
  METHOD <- paste(if(Wald) "Wald" else "Likelihood Ratio", "Test of association")
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, 
               conf.int = CINT, alternative = alternative,
               method = METHOD, 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

hazard.ratio.test(time = dat$time, status = dat$status, group = dat$receptor, rev = TRUE)
hazard.ratio.test(time = dat$time, status = dat$status, group = dat$receptor, Wald = FALSE, rev = TRUE)

# Fig 10.7
library(flexsurv)
fite <- flexsurvreg(formula = Surv(time, status) ~ receptor, data = breast.survival, dist="exp")
plot(fite, ci=FALSE)

# 10.2.2
# Exact conditional methods for a single 1 x 2 table

d1 <- sum(breast.survival$status[breast.survival$receptor=="Low" & breast.survival$stage=="III"])
m <- sum(breast.survival$status[breast.survival$stage=="III"])
n1 <- sum(breast.survival$time[breast.survival$receptor=="Low" & breast.survival$stage=="III"])
d2 <- m - d1
n <- sum(breast.survival$time[breast.survival$stage=="III"])
n2 <- n - n1



alpha <- (1-0.95)/2
1 - pbinom(12-1, size = 20, prob = 0.361)
pbinom(12, size = 20, prob = 0.809)


f1 <- function(x)1 - pbinom(d1-1, size = m, prob = x) - alpha
f2 <- function(x)pbinom(d1, size = m, prob = x) - alpha
f1.out <- uniroot(f1, interval = c(0,1))
f2.out <- uniroot(f2, interval = c(0,1))
CINTpi <- c(f1.out$root, f2.out$root)
CINT <- CINTpi*n2/((1 - CINTpi)*n1) # HR CI

# exact test of association
n1 <- sum(breast.survival$time[breast.survival$receptor=="Low" & breast.survival$stage=="III"])
n <- sum(breast.survival$time[breast.survival$stage=="III"])
pi_0 <- n1/n
# or
null <- 1
pi_0 <- (null * n1) / (null * n1 + n2)

min( min(pbinom(q = d1, size = m, prob = pi_0), 1 - pbinom(q = d1-1, size = m, prob = pi_0)) * 2, 1)

exact.ratio.test <- function(time, status, group, null=1, conf.level = 0.95, rev=FALSE){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  g <- levels(group)
  if(rev) g <- rev(g)
  d1 <- sum(status[group==g[1]])
  m <- sum(status)
  d2<- m - d1
  n1 <- sum(time[group==g[1]])
  n <- sum(time)
  n2 <- n - n1
  # estimate
  if(d1 == 0 || d2 == 0) est <- (d1 + 0.5) * n2 / (d2  + 0.5) * n1
  else est <- (d1 * n2) / (d2 * n1)
  names(est) <- "hazard ratio"
  names(null) <- names(est)
  # CI
  alpha <- (1-conf.level)/2
  f1 <- function(x)1 - pbinom(d1-1, size = m, prob = x) - alpha
  f2 <- function(x)pbinom(d1, size = m, prob = x) - alpha
  f1.out <- uniroot(f1, interval = c(0,1))
  f2.out <- uniroot(f2, interval = c(0,1))
  CINTpi <- c(f1.out$root, f2.out$root)
  CINT <- CINTpi*n2/((1 - CINTpi)*n1) # HR CI
  attr(CINT, "conf.level") <- conf.level
  
  # test
  # exact test of association
  pi_0 <- (null * n1) / (null * n1 + n2)
  p.value <- min( min(pbinom(q = d1, size = m, prob = pi_0), 1 - pbinom(q = d1-1, size = m, prob = pi_0)) * 2, 1)
  RVAL <- list(p.value = p.value, estimate = est, null.value = null,
               conf.int = CINT, alternative = alternative,
               method = "Exact Hazard Ratio Test of Association", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

exact.ratio.test(time = dat$time, status = dat$status, group = dat$receptor, rev = TRUE)

# 10.2.3
# Mantel-Haenszel test of association for person-time data


# Example 10.11
dat <- breast.survival[breast.survival$stage=="III",]
n1 <- sum(dat$time[dat$receptor=="Low"])
n2 <- sum(dat$time[dat$receptor=="High"])
d1 <- sum(dat$status[dat$receptor=="Low"])
d2 <- sum(dat$status[dat$receptor=="High"])
hr.hat <- (d1 * n2)/(d2 * n1)
m <- hr.hat*(d2/n2)*n1 + (d2/n2)*n2
e1 <- (n1*m)/(n1 + n2)
e2 <- (n2*m)/(n1 + n2)

X_pt <- (d1 - e1)^2 / e1 + (d2 - e2)^2 / e2
pchisq(X_pt, df = 1, lower.tail = FALSE)

# just for two groups
mantelhaen.pt.test <- function(time, status, group, rev=FALSE){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  null <- 1
  g <- levels(group)
  if(rev) g <- rev(g)
  n1 <- sum(time[group==g[1]])
  n2 <- sum(time[group==g[2]])
  d1 <- sum(status[group==g[1]])
  d2 <- sum(status[group==g[2]])
  est <- (d1 * n2)/(d2 * n1)
  names(est) <- "hazard ratio"
  names(null) <- names(est)
  m <- est*(d2/n2)*n1 + (d2/n2)*n2
  e1 <- (n1*m)/(n1 + n2)
  e2 <- (n2*m)/(n1 + n2)
  STATISTIC <- (d1 - e1)^2 / e1 + (d2 - e2)^2 / e2
  if(e1 < 5 || e2 < 5) warning(paste("expected counts less than 5: e1 =",e1,", e2 =",e2))
  names(STATISTIC) <- "Mantel-Haenszel X-squared"
  p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, null.value = null,
               alternative = alternative,
               method = "Mantel-Haenszel test of association for person-time data", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

mantelhaen.pt.test(time = dat$time, status = dat$status, group = dat$receptor, rev = TRUE)
# test warning
dat2 <- dat[-(1:8),]
mantelhaen.pt.test(time = dat2$time, status = dat2$status, group = dat2$receptor, rev = TRUE)

# Example 10.12
library(epitools)
rateratio(x = c(6,23,20), y = c(3637, 4797, 1037), method = "wald")
# disease counts
dc <- with(breast.survival, tapply(status, stage, sum))
# person time
ptime <- with(breast.survival, tapply(time, stage, sum))
rateratio(x = as.vector(dc), y = as.vector(ptime), method = "wald")

# mantel-haenszel test for 1 x I table
d <- with(breast.survival, tapply(status, stage, sum))
n <- with(breast.survival, tapply(time, stage, sum))
e <- n * sum(d) / sum(n)
X_pt <- sum((d - e)^2 / e)
pchisq(X_pt, df = length(d) - 1, lower.tail = FALSE)

# hazard ratios?
(d*n[1])/(d[1]*n)

# two or more groups
mantelhaen.pt.test <- function(time, status, group){
  dname <- deparse(substitute(time))
  alternative <- "two.sided"
  null <- 1
  d <- tapply(status, group, sum)
  n <- tapply(time, group, sum)
  est <- (d*n[1])/(d[1]*n)
  names(est) <- levels(group)
  names(null) <- "hazard ratio"
  e <- n * sum(d) / sum(n)
  STATISTIC <- sum((d - e)^2 / e)
  p.value <- pchisq(STATISTIC, df = length(d) - 1, lower.tail = FALSE)
  if(any(e < 5)) warning(paste("expected counts less than 5:",e[e < 5]))
  names(STATISTIC) <- "Mantel-Haenszel X-squared"
  p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  RVAL <- list(statistic = STATISTIC, p.value = p.value, estimate = est, null.value = null,
               alternative = alternative,
               method = "Mantel-Haenszel test of association for person-time data", 
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

mantelhaen.pt.test(time = breast.survival$time, status = breast.survival$status, group = breast.survival$stage)
with(breast.survival, mantelhaen.pt.test(time, status, stage))

mantelhaen.pt.test(time = dat$time, status = dat$status, group = dat$receptor)
# test warning
dat2 <- dat[-(1:8),]
mantelhaen.pt.test(time = dat2$time, status = dat2$status, group = dat2$receptor)
