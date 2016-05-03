# Biostats in epidemiology
# ch 4 work


# 12.2 --------------------------------------------------------------------

# Example 12.2

library(epitools)
schizophrenia

# crude death rate and directly standardized death rate
with(schizophrenia, 
     ageadjust.direct(count = cohort.deaths, pop = cohort.py, stdpop = alberta.pop)) * 10^3

# crude death rates
Ra <- sum((schizophrenia$cohort.py/sum(schizophrenia$cohort.py)) * 
      (schizophrenia$cohort.deaths/schizophrenia$cohort.py))
Rb <- sum((schizophrenia$alberta.pop/sum(schizophrenia$alberta.pop)) * 
      (schizophrenia$alberta.deaths/schizophrenia$alberta.pop))

# directly standardized death rate
Ras <- sum((schizophrenia$alberta.pop/sum(schizophrenia$alberta.pop)) * 
      (schizophrenia$cohort.deaths/schizophrenia$cohort.py))
Rbs <- sum((schizophrenia$alberta.pop/sum(schizophrenia$alberta.pop)) * 
      (schizophrenia$alberta.deaths/schizophrenia$alberta.pop))
# same as Rb

# crude rate ratio
CRR <- Ra/Rb

# standardized rate ratio
SRR <- Ras/Rb

# var of standardized death rate
Nsk <- schizophrenia$alberta.pop
Ns <- sum(schizophrenia$alberta.pop)
Dak <- schizophrenia$cohort.deaths
Nak <- schizophrenia$cohort.py

Dbk <- schizophrenia$alberta.deaths
Nbk <- schizophrenia$alberta.pop

# sqrt(sum((Nsk/Ns)^2 * (Dak/(Nak^2))))
varRas <- sum((Nsk/Ns)^2 * (Dak/(Nak^2)))
varRbs <- sum((Nsk/Ns)^2 * (Dbk/(Nbk^2)))
varlogSRR <- varRas/(Ras^2) + varRbs/(Rb^2)

# 95 CI for SRR
exp(log(SRR) + c(-1,1) * qnorm(0.975) * sqrt(varlogSRR))


# function for this
std.rate.ratio <- function(count, pop, stdcount, stdpop, conf.level=0.95){
  Ra <- sum((pop/sum(pop)) * 
              (count/pop))
  Rb <- sum((stdpop/sum(stdpop)) * 
              (stdcount/stdpop))
  
  # directly standardized death rate
  Ras <- sum((stdpop/sum(stdpop)) * 
               (count/pop))

  # crude rate ratio
  CRR <- Ra/Rb
  
  # standardized rate ratio
  SRR <- Ras/Rb
  
  # var of standardized death rate
  Nsk <- stdpop
  Ns <- sum(stdpop)
  Dak <- count
  Nak <- pop
  
  Dbk <- stdcount
  Nbk <- stdpop
  
  # sqrt(sum((Nsk/Ns)^2 * (Dak/(Nak^2))))
  varRas <- sum((Nsk/Ns)^2 * (Dak/(Nak^2)))
  varRbs <- sum((Nsk/Ns)^2 * (Dbk/(Nbk^2)))
  varlogSRR <- varRas/(Ras^2) + varRbs/(Rb^2)
  
  # 95 CI for SRR
  alpha <- (1 - conf.level)/2
  CINT <- exp(log(SRR) + c(-1,1) * qnorm(1 - alpha) * sqrt(varlogSRR))
  list("Standardize rate ratio" = SRR, "95% CI for SRR" = CINT,
       "Crude rate ratio" = CRR, 
       "directly standardized death rate" = Ras)
}

with(schizophrenia, std.rate.ratio(count = cohort.deaths, pop = cohort.py,
                                   stdcount = alberta.deaths, stdpop = alberta.pop))


# 12.3 --------------------------------------------------------------------

# Example 12.3
str(schizophrenia)
Da <- sum(schizophrenia$cohort.deaths)
Rsk <- schizophrenia$alberta.deaths/schizophrenia$alberta.pop
Nak <- schizophrenia$cohort.py
Ea <- sum(Rsk*Nak)
# standardized mortality ratio
SMR <- Da/Ea

varSMR <- SMR/Ea

# 95% CI
SMR + c(-1,1)*qnorm(0.975)*sqrt(varSMR)

# This returns a slightly different CI
ageadjust.indirect(count = schizophrenia$cohort.deaths, 
                   pop = schizophrenia$cohort.py, 
                   stdcount = schizophrenia$alberta.deaths, 
                   stdpop = schizophrenia$alberta.pop)

# exact CI (if Ea < 5)
# use 10.9 and 10.10 with n = Ea and d = Da (time = Ea and status = Da)
exact.rate.test(time = Ea, status = Da)$conf.int
exact.rate.test(time = Ea, status = Da)
# For the 10-19 age group with Da = 2 and Ea = 0.376
exact.rate.test(time = 0.376, status = 2)

# hypothesis of no mortality difference between cohort and the standard pop'n:
STATISTIC <- ((Da - Ea)^2)/Ea


smr.test <- function(count, pop, stdcount, stdpop, conf.level = 0.95){
  dname <- paste("\n  ", deparse(substitute(count)), 
                 "\n  ", deparse(substitute(pop)),
                 "\n  ", deparse(substitute(stdcount)),
                 "\n  ", deparse(substitute(stdpop)))
  alternative <- "two.sided"
  
  Da <- sum(count)
  Rsk <- stdcount/stdpop
  Nak <- pop
  Ea <- sum(Rsk*Nak)
  # standardized mortality ratio
  est <- Da/Ea
  names(est) <- "Standardized Mortality Ratio"
  null <- 1
  names(null) <- names(est)

  if(Ea < 5){
    CINT <- exact.rate.test(time = Ea, status = Da)$conf.int
    attr(CINT, "conf.level") <- conf.level
    p.value <- exact.rate.test(time = Ea, status = Da)$p.value
    METHOD <- paste("Exact test of no mortality difference between cohort and standard population")
    RVAL <- list(p.value = p.value, 
                 estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = METHOD, 
                 data.name = dname) 
    
  } else {
    varSMR <- est/Ea
    # 95% CI
    alpha <- (1-conf.level)/2
    CINT <- est + c(-1,1)*qnorm(1 - alpha)*sqrt(varSMR)
    attr(CINT, "conf.level") <- conf.level
    STATISTIC <- ((Da - Ea)^2)/Ea
    p.value <- pchisq(q = STATISTIC, df = 1, lower.tail = FALSE)
    names(STATISTIC) <- "X-squared"
    METHOD <- paste("Test of no mortality difference between cohort and standard population")
    RVAL <- list(statistic = STATISTIC, parameter = c(df = 1), p.value = p.value, 
                 estimate = est, null.value = null,
                 conf.int = CINT, alternative = alternative,
                 method = METHOD, 
                 data.name = dname) 
    
  }
  class(RVAL) <- "htest"
  return(RVAL)
}

with(schizophrenia, 
     smr.test(count = cohort.deaths,
                    pop = cohort.py, 
                    stdcount = alberta.deaths, 
                    stdpop = alberta.pop))

# Exact test (expected number of deaths < 5)
with(subset(schizophrenia,age.group=="10-19"), 
     smr.test(count = cohort.deaths,
                    pop = cohort.py, 
                    stdcount = alberta.deaths, 
                    stdpop = alberta.pop))

# top <- schizophrenia[,1:3] 
# top$exp <- 1
# bot <- schizophrenia[,c(1,4,5)]
# bot$exp <- 2
# names(top) <- names(bot) <- c("age.group","deaths","pop", "exp")
# dat <- rbind(top,bot)
# 
# schizL <- melt(schizophrenia, id.vars = c("age.group","cohort.deaths"))
# 
# lrt.homogeneity(time, status, exposure, strata)


# 12.4 --------------------------------------------------------------------


