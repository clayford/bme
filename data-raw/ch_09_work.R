# Biostats in epidemiology
# ch 9 work

library(survival)

# Table 9.1
# ILow <- c(50,51,51,53,53,54,54,55,56,56,57,60)
# ILowc <- c(0,1,0,0,0,0,0,0,1,0,0,0)
# 
# IHigh <- c(10,34,34,47,47,49,49,rep(50,7),rep(51,6), rep(52,5), rep(53,6), rep(54,5),
#            55,55,56,56,rep(57,5), rep(58,5), rep(59,4), 60, 60, 60)
# IHighc <- c(1,1,0,1,1,rep(0, sum(c(2,7,6,5,6,5,2,2,5,5,4,3))))
# 
# IILow <- c(4,9,13,21,29,29,40,46,49,49,52,52,53,54,55,55,56,57,57,58,58,59,60)
# IILowc <- c(0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
#             1,0,0,0,0,0)
# 
# IIHigh <- c(11,16,21,23,23,24,33,33,36,36,36,37,45,46,49,49,rep(50,6),
#             rep(51,4),rep(52,5), rep(53,5), rep(54,4), rep(55,4), rep(56,6), rep(57,4), 58, 58,
#             rep(58,8), rep(59,5), rep(60,6))
# IIHighc <- c(rep(1,10),0,1,1,1,rep(0, sum(c(2,6,4,5,5,4,4,6,4))), 1, 1,
#              rep(0, sum(c(8,5,6))))
# 
# IIILow <- c(9, 12, 14, 15, 15, 17, 21, 22, 23, 23, 31, 34, 35, 53, 60)
# IIILowc <- c(rep(1,4), 0, rep(1,8), 0, 0)
# 
# IIIHigh <- c(7, 9, 17, 21, 22, 22, 34, 34, 41, 49, 52, 55, 56, 58, 58, 59, 59)
# IIIHighc <- c(0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, rep(0,5))
# 
# r0 <- c(12, 57, 23, 75, 15, 17)
# tab9.1 <- data.frame(time = c(ILow, IHigh, IILow, IIHigh, IIILow, IIIHigh),
#                      status = c(ILowc, IHighc, IILowc, IIHighc, IIILowc, IIIHighc),
#                      stage = c(rep("I", (12 + 57)), rep("II", (23 + 75)), rep("III", (15 + 17))),
#                      receptor = rep(rep(c("Low","High"), 3), r0))
# tab9.1 <- tab9.1[order(tab9.1$time),]
# 
# write.csv(tab9.1, file = "BreastCancer.csv", row.names = FALSE)

breast.survival <- read.csv("BreastCancer.csv")


# 9.1 ---------------------------------------------------------------------

breast.survival$time

# death times (tau_j)
tau_j <- unique(breast.survival$time[breast.survival$status == 1])

breast.survival$intervals <- cut(breast.survival$time, breaks = c(0, tau_j, Inf), right = FALSE)

# the risk set: the r_j column in Table 9.2
sum(breast.survival$time >= tau_j[1])
sum(breast.survival$time >= tau_j[2])
sum(breast.survival$time >= tau_j[3])
sum(breast.survival$time >= tau_j[4])

# get all risk sets (r_j)
r_j <- sapply(tau_j, function(x)sum(breast.survival$time >= x))

# number of deaths at tau_j (a_j)
sum(breast.survival$time[breast.survival$status == 1] == tau_j[1])
sum(breast.survival$time[breast.survival$status == 1] == tau_j[2])
sum(breast.survival$time[breast.survival$status == 1] == tau_j[3])
sum(breast.survival$time[breast.survival$status == 1] == tau_j[4])
sum(breast.survival$time[breast.survival$status == 1] == tau_j[9])
sum(breast.survival$time[breast.survival$status == 1] == tau_j[10])
sum(breast.survival$time[breast.survival$status == 1] == tau_j[11])
sum(breast.survival$time[breast.survival$status == 1] == tau_j[12])

# get all number of deaths at tau_j
a_j <- sapply(tau_j, function(x)sum(breast.survival$time[breast.survival$status==1] == x))

# conditional probability of surviving to just after tau_j (p_j)
p_j <- (r_j - a_j)/r_j

#probability of surviving from tau_0 to tau_j, the Kaplan-Meier estimate (S_j)
S_j <- cumprod(p_j)

# estimate of var(S_j): greenwood
sqrt((S_j^2) * cumsum((1 - p_j)/(p_j*r_j)))

# estimate of var(S_j): Kalbfleisch-Prentice method
(1/(log(S_j)^2))*cumsum((1 - p_j)/(p_j*r_j))



# top of page 177 work for Kalbfleisch-Prentice method
sqrt((1 - p_j[1])/(p_j[1]*r_j[1]) + (1 - p_j[2])/(p_j[2]*r_j[2]))

sqrt(((1 - p_j[1])/(p_j[1]*r_j[1]) + (1 - p_j[2])/(p_j[2]*r_j[2])) / log(S_j[2])^2)

Su <- log(-log(S_j[2])) - sqrt((((1 - p_j[1])/(p_j[1]*r_j[1]) + (1 - p_j[2])/(p_j[2]*r_j[2])) / log(S_j[2])^2)) * qnorm(0.975)
exp(-exp(Su))


# Table 9.2 (without CIs) 
tab9.2 <- data.frame(j = 1:30, tau_j, a_j, r_j, p_j, S_j)

# Example 9.1

# create a Surv object
SurvBC <- Surv(breast.survival$time, breast.survival$status)
surv.all <- survfit(SurvBC ~ 1, conf.type = "log-log")
# Table 9.2 again with CIs
summary(surv.all)
# Note: the standard error column is based on greenwood's formula and not used
# in the "log-log" CIs
sqrt((S_j^2) * cumsum((1 - p_j)/(p_j*r_j)))



# Figure 9.2
plot(surv.all, ylim=c(0.5,1.0), ylab = "Survival probability", xlab = "Time")


# 9.2 ---------------------------------------------------------------------



# Example 9.2
surv.all2 <- survfit(SurvBC ~ receptor, data=breast.survival)
# Figure 9.3(a)
plot(surv.all2, lty = 2:3, ylim=c(0.4,1.0), ylab = "Survival probability", xlab = "Time", 
     axes = FALSE)
axis(side = 1, at = seq(0,60,12))
axis(side = 2, at = seq(0.4,1.0,0.2))
text(42,0.9,"High level")
text(42,0.7,"Low level")

# logrank test
sdout <- survdiff(SurvBC ~ receptor, data=breast.survival)
sdout

sum.out <- summary(surv.all2)
str(sum.out)

# Table 9.4
# death times (tau_j)
tau_j <- unique(breast.survival$time[breast.survival$status == 1])

# get all number of deaths at tau_j for recepter level low (a_1j)
a_1j <- sapply(tau_j, function(x)sum(breast.survival$time[breast.survival$status==1 & breast.survival$receptor == "Low"] == x))
# get all number of deaths at tau_j for recepter level high (a_2j)
a_2j <- sapply(tau_j, function(x)sum(breast.survival$time[breast.survival$status==1 & breast.survival$receptor == "High"] == x))

# get all risk sets (r_1j) for receptor level low
r_1j <- sapply(tau_j, function(x)sum(breast.survival$time[breast.survival$receptor == "Low"] >= x))
r_2j <- sapply(tau_j, function(x)sum(breast.survival$time[breast.survival$receptor == "High"] >= x))

# e.hat_1j
e.hat_1j <- (r_1j*(a_1j + a_2j))/(r_1j + r_2j)

# v.hat_0j
v.hat_0j <- (r_1j * r_2j * (a_1j + a_2j) * ((r_1j - a_1j) + (r_2j - a_2j))) / ((r_1j + r_2j)^2 * ((r_1j + r_2j) - 1))

tab9.4 <- data.frame(j = 1:30, tau_j, a_1j, a_2j, r_1j, r_2j, e.hat_1j, v.hat_0j)

# logrank test statistic using Table 9.4
X_mh <- ((sum(a_1j) - sum(e.hat_1j))^2) / sum(v.hat_0j)
pchisq(X_mh, df = 1, lower.tail = FALSE)

or_analysis <- data.frame(a_1j, a_2j, 
                          b_1j = r_1j - a_1j, b_2j = r_2j - a_2j, 
                          tau_j)
library(reshape2)
or_analysisL <- melt(or_analysis, id.vars = "tau_j", value.name = "Freq")
library(tidyr)
or_analysisL <- separate(data = or_analysisL, col = variable, into = c("Survival","Receptor"))
or_analysisL$Survival <- factor(or_analysisL$Survival, labels = c("dead","alive"))
or_analysisL$Receptor <- factor(or_analysisL$Receptor, labels = c("low","high"))
or_analysisT <- xtabs(Freq ~ Survival + Receptor + tau_j, data = or_analysisL)
mantelhaen.test(or_analysisT)

# Function for odds ratio

fo <- as.formula("Freq ~ Survival + Receptor + tau_j")

xtabs(fo, data = or_analysisL)

xtabs2 <- function(formula, data){
  xtabs(formula, data = data)
}

xtabs2(formula = Freq ~ Survival + Receptor + tau_j, data = or_analysisL)

# function: odds ratio for censored survival data
or.csd <- function(data, time, censor, var, rev=FALSE){
  if(!rev)levs <- levels(data[,var])
  else levs <- levels(data[,var])[c(2,1)]
  # death times (tau_j)
  tau_j <-   unique(data[,time][data[,censor] == 1])
  
  # get all number of deaths at tau_j for level 1
  a_1j <- sapply(tau_j, function(x)sum(data[,time][data[,censor]==1 & data[,var] == levs[1]] == x))
  # get all number of deaths at tau_j for level 2
  a_2j <- sapply(tau_j, function(x)sum(data[,time][data[,censor]==1 & data[,var] == levs[2]] == x))
  
  # get all risk sets (r_1j) 
  r_1j <- sapply(tau_j, function(x)sum(data[,time][data[,var] == levs[1]] >= x))
  r_2j <- sapply(tau_j, function(x)sum(data[,time][data[,var] == levs[2]] >= x))
  
  or_analysis <- data.frame(a_1j, a_2j, 
                            b_1j = r_1j - a_1j, b_2j = r_2j - a_2j, 
                            tau_j)
  or_analysisL <- reshape2::melt(or_analysis, id.vars = "tau_j", value.name = "Freq")
  or_analysisL <- tidyr::separate(data = or_analysisL, col = variable, into = c("Survival",var))
  or_analysisL$Survival <- factor(or_analysisL$Survival, labels = c("dead","alive"))
  or_analysisL[,var] <- factor(or_analysisL[,var], labels = c(levs[1],levs[2]))
  fo <- as.formula(paste("Freq ~ Survival +", var,"+ tau_j"))
  or_analysisT <- xtabs(fo, data = or_analysisL)
  mantelhaen.test(or_analysisT)
}
or.csd(data = breast.survival, time = "time", censor = "status", var = "receptor", rev = TRUE)



test1 <- function(data, time, censor, var){
  # death times (tau_j)
  unique(data[,time][data[,censor] == 1])
  
}

test1(data = breast.survival, time = "time", censor = "status", var = "receptor")
  
# 9.2.2
# assessment of proportional hazards assumption
# Figure 9.3(b)
sum.out <- summary(surv.all2)
plot(sum.out$time[sum.out$strata=="receptor=Low"], 
     log(-log(sum.out$surv[sum.out$strata == "receptor=Low"])), 
     type = "s", ylim=c(-5, 0), xlim=c(0,60), xlab = "Time", ylab = "Log - minus - Log")
lines(sum.out$time[sum.out$strata=="receptor=High"], 
     log(-log(sum.out$surv[sum.out$strata == "receptor=High"])), 
     type = "s", lty=2)

# test for linear trend
prop.trend.test(x = a_1j, n = r_1j)

# using Cox regression (adapting example from documentation)
fit <- coxph(Surv(time, status) ~ receptor,  
             data=breast.survival) 
temp <- cox.zph(fit) 
print(temp)                  # display the results 
plot(temp)                   # plot curves 


# Table 9.6
# data from a cohort study of women with stage II or stage IIIA ovarian cancer.
# grade is an indicator of the malignant potential of the tumor.
ovarian <- data.frame(grade = rep(c("Low","High"), c(15,20)), 
                     time = c(28, 89, 175, 195, 309, 377, 393, 421, 447, 462, 709, 744, 770,
                              1106, 1206,
                              34, 88, 137, 199, 280, 291, 299, 300, 309, 351, 358, 369, 369, 370,
                              375, 382, 392, 429, 451, 1119),
                     status = c(rep(1,5), rep(0,4), 1, rep(0,5),
                                rep(1, 6), 0, 0, rep(1,9), 0, 1, 0))
write.csv(ovarian, file = "ovarian.csv", row.names = FALSE)
# Example 9.5

fit <- survfit(Surv(time, status) ~ grade, data = ovarian) 
# Fig 9.4(a)
plot(fit, lty = 2:3, xlim=c(0,500)) 

# logrank test
survdiff(Surv(time, status) ~ grade, data=ovarian)
# evidence of mortality difference between two groups

# assess proportional hazard assumption
# Figure 9.4(b)
sum.out <- summary(fit)
ylim <- round(range(log(-log(sum.out$surv))))
xlim <- c(0, 500)
plot(sum.out$time[sum.out$strata=="grade=Low"], 
     log(-log(sum.out$surv[sum.out$strata == "grade=Low"])), 
     type = "s", ylim=ylim, xlim = xlim, xlab = "Time", ylab = "Log - minus - Log")
lines(sum.out$time[sum.out$strata=="grade=High"], 
      log(-log(sum.out$surv[sum.out$strata == "grade=High"])), 
      type = "s", lty=2)

or.csd(data = ovarian, time = "time", censor = "status", var = "grade")

# 9.2.3
# Example 9.6
fit <- survfit(Surv(time, status) ~ stage, data=breast.survival)
# fig 9.5(a)
plot(fit, lty = 1:3, ylim = c(0.2, 1.0))
# log rank test of association
survdiff(Surv(time, status) ~ stage, data=breast.survival)

# log-minus-log curve
# fig 9.5(b)
sum.out <- summary(fit)
ylim <- c(-5,1)
xlim <- c(0, 60)
plot(sum.out$time[sum.out$strata=="stage=I"], 
     log(-log(sum.out$surv[sum.out$strata == "stage=I"])), 
     type = "s", ylim=ylim, xlim = xlim, xlab = "Time", ylab = "Log - minus - Log")
lines(sum.out$time[sum.out$strata=="stage=II"], 
      log(-log(sum.out$surv[sum.out$strata == "stage=II"])), 
      type = "s", lty=2)
lines(sum.out$time[sum.out$strata=="stage=III"], 
      log(-log(sum.out$surv[sum.out$strata == "stage=III"])), 
      type = "s", lty=3)
legend("topleft", legend = c("Stage I", "Stage II", "Stage III"), lty = 1:3)

# Fig 9.6
breast.survival$stg.recep <- with(breast.survival, interaction(stage, receptor))
fit <- survfit(Surv(time, status) ~ stg.recep, data=breast.survival)
plot(fit, lty = 1:6, xlab = "Time", ylab = "Survival probability")
legend("bottomleft", legend = levels(breast.survival$stg.recep), lty = 1:6)



# 9.3 Actuarial method

# create table 9.13
# death times (tau_j) intervals
tau_j <- seq(0,60,12)

# add intervals to tab9.1
breast.survival$interval <- cut(breast.survival$time, breaks = c(0,seq(12,60,12),Inf), right = FALSE)

# get all risk sets (r_j) for each interval
r_j <- sapply(tau_j, function(x)sum(breast.survival$time >= x))

# get all number of deaths at tau_j intervals
a_j <- with(breast.survival, tapply(status, interval, sum))

# get all number of censored at tau_j intervals
c_j <- with(breast.survival, tapply(status, interval, function(x)(sum(x == 0))))

# r`_j: the effective sample size
rp_j <- r_j - (c_j/2)

# conditional probability of surviving 
p_j <- 1 - (a_j/rp_j)

#probability of surviving from tau_0 to tau_j, the Kaplan-Meier estimate (S_j)
S_j <- cumprod(p_j)

tab9.13 <- data.frame(tau_j, a_j, r_j, c_j, rp_j, p_j, S_j = c(1.0, S_j[-6]))

plot(tab9.13$tau_j, tab9.13$S_j, type = "b", xlab = "Time", ylab = "Survival probability", 
     axes = FALSE, ylim = c(0.6, 1.0))
axis(side = 1, at = seq(0,60,12))
axis(side = 2, at = seq(0.6,1.0,0.1))
