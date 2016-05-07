# Biostats in epidemiology
# ch 14 work


# 14.2 --------------------------------------------------------------------

# sample size for closed cohort study
ss.closed.chort <- function(pi2, delta, sig.level = 0.05, power = 0.80, rho = 1, 
                            type=c("RD","RR","OR")){
  type <- match.arg(type)
  alternative <- "two.sided"
  if(type == "RD") pi1 <- pi2 + delta
  else if(type == "RR") pi1 <- pi2*delta
  else pi1 <- delta*pi2 / (delta*pi2 + (1 - pi2))
  pi0 <- (pi1 + (pi2 * rho)) / (1 + rho)
  num <- qnorm(1 - sig.level/2) * sqrt(pi0 * (1 - pi0) * ((1 + rho)/ rho)) + 
  qnorm(1 - (1 - power)) * sqrt(pi1 * (1 - pi1) + (pi2 * (1 - pi2))/rho)
  r1 <- ceiling((num^2) / (delta^2))
  r2 <- r1 * rho
  NOTE <- "r1 is the number of exposed subjects needed for the study"
  METHOD <- paste(switch(type, RD = "Risk Difference", RR = "Risk Ratio", OR = "Odds Ratio"), 
                         "- Sample size for closed cohort study")
  structure(list(r1 = r1, r2 = r2, delta = delta, rho = rho, sig.level = sig.level, 
                 power = power, alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}

# Example 14.2
ss.closed.chort(pi2 = 0.05, delta = c(0.01, 0.05, 0.10, 0.20, 0.30), power = 0.8, rho = 1)

# Example 14.3
ss.closed.chort(pi2 = 0.05, delta = 0.05, power = 0.8, rho = c(1:5,10,20))

# Risk Difference
ss.closed.chort(pi2 = 0.05, delta = 0.05, power = 0.8, rho = 1)
# Risk Ratio
ss.closed.chort(pi2 = 0.05, delta = 0.05, power = 0.8, rho = 1, type = "RR")
# Odds Ratio
ss.closed.chort(pi2 = 0.05, delta = 0.05, power = 0.8, rho = 1, type = "OR")



# 14.6 --------------------------------------------------------------------

# Power

rho <- 2
p2 <- 0.10
OR <- 2
m <- 100 # cases
p1 <- (OR * p2)/(OR*p2 + (1 - p2))
p0 <- (p1 + p2*rho)/(1 + rho)
p.body <- ((sqrt(m) * abs(p1 - p2)) - qnorm(0.975)*sqrt(p0*(1 - p0)*((1 + rho)/rho))) /
  sqrt(p1*(1 - p1) + ((p2*(1 - p2))/rho))
pnorm(p.body) # power

p.body <- quote({
  p1 <- (OR * p2)/(OR*p2 + (1 - p2))
  p0 <- (p1 + p2*rho)/(1 + rho)
  pnorm(((sqrt(m) * abs(p1 - p2)) - qnorm(0.975)*sqrt(p0*(1 - p0)*((1 + rho)/rho))) /
    sqrt(p1*(1 - p1) + ((p2*(1 - p2))/rho)))
  })
eval(p.body)

# find n
power <- 0.8
m <- NULL
uniroot(function(m) eval(p.body) - power, c(2, 1e+07), 
             tol = .Machine$double.eps^0.25, extendInt = "upX")$root


# Unmatched case control study
# Using power.t.test as a template
# m1 = number of cases
# m2 = rho * m1 = number of controls
# p1 = probability that a case has a history of exposure
# p2 = probability that a control has a history of exposure

# Odds Ratio for an incidence case-control study is the same whether we consider
# the row or column marginal totals fixes.

power.or.ucc <- function (n = NULL, power = NULL, rho = 1, p2, OR, sig.level = 0.05,  
          tol = .Machine$double.eps^0.25) 
{
  if (sum(sapply(list(n, power), is.null)) != 
      1) 
    stop("exactly one of 'n', and 'power'must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
    stop("'sig.level' must be numeric in [0, 1]")
  alternative <- "two.sided"
  p1 <- (OR * p2)/(OR*p2 + (1 - p2))
  p0 <- (p1 + p2*rho)/(1 + rho)
  p.body <- quote({
    p1 <- (OR * p2)/(OR*p2 + (1 - p2))
    p0 <- (p1 + p2*rho)/(1 + rho)
    pnorm(((sqrt(n) * abs(p1 - p2)) - qnorm(1 - (sig.level/2))*sqrt(p0*(1 - p0)*((1 + rho)/rho))) /
            sqrt(p1*(1 - p1) + ((p2*(1 - p2))/rho)))
  })
  
  if (is.null(power)) 
    power <- eval(p.body)
  else if (is.null(n)) 
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07), 
                 tol = .Machine$double.eps^0.25, extendInt = "upX")$root
  else stop("internal error", domain = NA)
  NOTE <- "n is number of cases; multiply by rho to obtain number of controls"
  METHOD <- "Unmatched case-control power calculation"
  structure(list(n = n, power = power, p2 = p2, OR = OR, rho = rho, sig.level = sig.level, 
                  alternative = alternative, note = NOTE, 
                 method = METHOD), class = "power.htest")
}


# Example 14.5 (Table 14.4)
sapply(c(2,3,4,5,10), function(x)power.icc(power = 0.8, rho = 1, p2 = 0.05, OR = x)$n)

# Example 14.7
power.icc(n = 100, rho = 2, OR = 2, p2 = 0.10)

# Example 14.8
power.icc(power = 0.8, rho = 6.73, OR = 3, p2 = 0.048)
