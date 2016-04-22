# Biostats in epidemiology
# ch 4 work



# Section 4.1 -------------------------------------------------------------

# log-odds trransformation

a1 <- rbinom(1e5, size = 10, prob = 0.4)
a2 <- rbinom(1e5, size = 25, prob = 0.2)
or <- (a1*(25 - a2)) / (a2*(10 - a1))
log.or <- log((a1*(25 - a2))) - log((a2*(10 - a1)))

# figs 4.1(a) and 4.1(b)
hist(or[or < 12.25], freq = FALSE, breaks = 10, main="Fig 4.1(a)")
br <- quantile(log.or, probs = c(0.01, 0.99), na.rm = TRUE)
hist(log.or[log.or > br[1] & log.or < br[2]], freq = FALSE, breaks = 10, main="Fig4.1(b)")


a1 <- rbinom(1e5, size = 25, prob = 0.4)
a2 <- rbinom(1e5, size = 50, prob = 0.2)
or <- (a1*(50 - a2)) / (a2*(25 - a1))
log.or <- log((a1*(25 - a2))) - log((a2*(10 - a1)))

# figs 4.2(a) and 4.2(b)
hist(or[or < 12.25], freq = FALSE, breaks = 10, main="Fig 4.2(a)")
br <- quantile(log.or, probs = c(0.01, 0.99), na.rm = TRUE)
hist(log.or[log.or > br[1] & log.or < br[2]], freq = FALSE, breaks = 10, main="Fig4.2(b)")

# log(OR) should be reasonably symmetric provided the means of the component
# distributions are 5 or more.

# confidence intervals

# example 4.1
# table 4.3
# antibody-diarrhea
antibody <- matrix(c(12,2,7,9),ncol=2, dimnames = list(diarrhea = c("yes","no"), antibody=c("low","high")))
antibody
addmargins(antibody)

binom.test(12, 14)
binom.test(7, 16)

library(epitools)
or.out <- oddsratio(t(antibody), rev="both", method="wald", verbose = T)
or.out

# table 4.4
expected(antibody)

# Pearson (provided in output of oddsratio as "chi.square")
chi_p <- sum(((antibody - expected(antibody))^2)/expected(antibody))
pchisq(chi_p, df = 1, lower.tail = FALSE)
# or do this
chisq.test(antibody, correct = F)

# Wald
chi_w <- (log(or.out$measure[2,1])^2)*(1/(sum(1/expected(antibody))))
pchisq(chi_w, df = 1, lower.tail = FALSE)

# likelihood ratio
chi_lr <- 2*sum(antibody*log(antibody/expected(antibody)))
pchisq(chi_lr, df = 1, lower.tail = FALSE)

# example 4.2
breast <- matrix(c(23,25,31,113),ncol=2, dimnames = list(survival = c("dead","alive"), 
                                                      receptor.level=c("low","high")))
breast
oddsratio(t(breast), re="both", method = "wald", verbose=TRUE)

# potential impact of misclassification
# solving equations at bottom of p. 99; displayed in Table 4.5(b)
lhs <- matrix(c(0.9, 0.01, 0.1, 0.99), nrow = 2, byrow = TRUE)
rhs <- breast[,1]
low <- solve(lhs, rhs)
rhs <- breast[,2]
high <- solve(lhs, rhs)
(low[1]*high[2])/(low[2]*high[1])

oddsratio(t(breastb), rev = "both", method = "wald")

oddsratio(breastb, method = "wald")$measure[2,1]
# function for missclassification

adjusted.odds.ratio <- function(data, sensitivity, specificity){
  if(!is.matrix(data) || nrow(data)!=2 || ncol(data) != 2) stop("data must be a 2 x 2 matrix")
  if(sensitivity > 1 || sensitivity < 0 || specificity > 1 || specificity < 0) stop("sensitivity and specificity must be in (0,1)")
  lhs <- matrix(c(sensitivity, 1-specificity, 1-sensitivity, specificity), nrow = 2, byrow = TRUE)
  rhs <- data[,1]
  low <- solve(lhs, rhs)
  rhs <- data[,2]
  high <- solve(lhs, rhs)
  adj <- (low[1]*high[2])/(low[2]*high[1])
  unadj <- (data[1,1]*data[2,2])/(data[1,2]*data[2,1])
  list("Adjusted odds ratio after accounting for misclassification"=adj,
       "Estimated odds ratio based on (possibly) misclassified data"= unadj)
}

adjusted.odds.ratio(data = breast, sensitivity = 0.9, specificity = 0.99)
adjusted.odds.ratio(data = breast.receptor, sensitivity = 0.9, specificity = 0.99)
adjusted.odds.ratio(data = breast, sensitivity = 0.9, specificity = 2)

# section 4.2 -------------------------------------------------------------

# exact conditional methods for a single 2 x 2 table

# 

# table 4.9
dhyper(x = 1:3, m = 3, n = 2, k = 3)

# example 4.3
t48 <- matrix(c(2,1,1,1), ncol = 2, dimnames = list(disease = c("yes", "no"), exposure = c("yes","no")))
addmargins(t48)
oddsratio(t(t48), rev = "both", method = "fisher")

# table 4.10
antibody
round(dhyper(x = 3:14, m = 19, n = 11, k = 14),4)*100
round(phyper(q = 3:14, m = 19, n = 11, k = 14),4)*100
round(phyper(q = (3:14) - 1, m = 19, n = 11, k = 14, lower.tail = FALSE),4)*100

# example 4.4
oddsratio(t(antibody), rev="both", method = "fisher")
fisher.test(antibody)
# getting the p-value "manually"
phyper(q = 11, m = 19, n = 11, k = 14, lower.tail = FALSE) + phyper(q = 5, m = 19, n = 11, k = 14)


# exact confidence intervals tend to be wider than asymptotic ones. 
# exact p-values are generally larger than their asymptotic counterparts.

# example 4.5
oddsratio(t(breast), rev = "both", method = "fisher")


# 4.3 ---------------------------------------------------------------------

# asymptotic conditional methods for a single 2 x 2 table

# example 4.6
oddsratio(t(t48), rev = "both", method = "fisher")

# example 4.7
oddsratio(t(antibody), rev = "both", method = "fisher")

MH <- ((12 - expected(antibody)[1])^2)/((14*16*11*19)/(30^2*29))
pchisq(MH, df = 1, lower.tail = FALSE)

# example 4.8
addmargins(breast)
oddsratio(t(breast), rev = "both", method = "fisher")
# explicit CI
r1 <- 48
r2<- 144
m1 <- 54
m2 <- 138
l <- max(0, r1 - m2) 
u <- min(r1, m1)
a1 <- 23
chat <- sum(choose(n = r1, k = l:u)*choose(r2, m1 - l:u)*3.329255^(l:u))
vhat <- sqrt(sum((l:u - a1)^2*choose(n = r1, k = l:u)*choose(r2, m1 - l:u)*3.329255^(l:u))/chat)
exp(log(3.329255) + c(-1,1)*qnorm(p = 0.975)/vhat)


# 4.6 ---------------------------------------------------------------------


choose(10,2)
1/20*45
t416 <- matrix(c(7,60,26,70,21,8),ncol=3,
              dimnames = list(survival=c("dead","alive"), stage=c("I","II","III")))
# table 4.16
addmargins(t416)
table.margins(t416)

# table 4.17
oddsratio(x = t(t416), rev = "columns", method = "wald")

library(coin)
cmh_test(as.table(t416))
chisq_test(as.table(t416))
chisq.test(x = t416) # same thing

ct <- chisq_test(as.table(t416),
                  scores = list("stage" = c(1,2,3)))
statistic(ct)^2


# table 4.18
t416e <- expected(t416)

# test for linear trend
# example 4.11 work
num <- sum(c(1,2,3)*(t416[1,] - t416e[1,]))^2
den <- sum(c(1,2,3)^2*t416e[1,]) - (sum(c(1,2,3)*t416e[1,])^2)/sum(t416[1,])
((sum(t416) - 1)/sum(t416[2,]))*(num/den)

# using the lbl_test in coin package
lbl_test(as.table(t416), scores = list("stage" = c(1,2,3)))
tout <- lbl_test(as.table(t416), scores = list("stage" = c(1,2,3)))
statistic(tout)^2

tout <- lbl_test(as.table(t416), scores = list("stage" = c(1,2,5)))
statistic(tout)^2

# see also the prop.trend.test in base R
prop.trend.test(x = c(7,26,21), n = c(67, 96, 29))
# slightly different answer