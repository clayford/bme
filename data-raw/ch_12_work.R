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
