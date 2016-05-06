# Biostats in epidemiology
# ch 13 work


# 13.1 --------------------------------------------------------------------

# Ordinary Life Table

# Example 13.1

Rj <- males$deaths.all/males$pop

# qj
1*Rj[1] / (1 + (1*Rj[1])/2)
4*Rj[2] / (1 + (4*Rj[2])/2)
5*Rj[3] / (1 + (5*Rj[3])/2)
# 85 (last age group)
5*Rj[19] / (1 + (5*Rj[19])/2)


nj <- diff(age.group)
qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
pj <- 1 - qj

# l(0) usually defined to be something like 100,000
l0 <- 1e5
lxj <- cumprod(c(l0,pj[-length(pj)]))

dj <- qj * lxj

Lj <- dj/Rj


Txj <- sapply(1:length(Lj), function(x)sum(Lj[x:length(Lj)]))

exj <- Txj/lxj

data.frame(xj = males$age.group, qj = round(qj,5), pj = round(pj, 5), 
           lxj = round(lxj), dj = round(dj), Lj = round(Lj), 
           Txj = round(Txj), exj = round(exj,2))


# function to create ordinary life table

olt <- function(group, deaths, pop, radix=1e5, phi=1){
  # annual death rate in population
  Rj <- (deaths/pop)
  # multiply by phi for sensitivity analysis
  Rj <- Rj * phi
  # length of age group
  nj <- diff(group)
  # conditional probability of dying 
  qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
  # conditional probability of surviving
  pj <- 1 - qj
  # the radix: number of individuals in the OLT birth cohort
  l0 <- radix
  # number of survivors to age x
  lxj <- cumprod(c(l0,pj[-length(pj)]))
  # number of deaths in age group
  dj <- qj * lxj
  # number of person years in age group
  Lj <- dj/Rj
  # Number of person years after age x
  Txj <- sapply(1:length(Lj), function(x)sum(Lj[x:length(Lj)]))
  # Life expectancy at age x
  exj <- Txj/lxj
  # collect results into data frame and display
  data.frame(group = group, qj = round(qj,5), pj = round(pj, 5), 
             lxj = round(lxj), dj = round(dj), Lj = round(Lj), 
             Txj = round(Txj), exj = round(exj,2))
}
with(males, olt(group = age.group, deaths = deaths.all, pop = pop))
with(males, olt(group = age.group, deaths = deaths.all, pop = pop, phi = 0.9))


# 13.2 --------------------------------------------------------------------

# multiple decrement life table
# makes it possible to examine specific causes of death
# example of a competing risks model

xj <- males$age.group

# step 1
Rj <- males$deaths.all/males$pop
nj <- diff(males$age.group)
qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
Djk <- males$neoplasm.deaths 
Dj <- males$deaths.all
qjk <- (Djk/Dj) * qj

# step 2
pj <- 1 - qj
l0 <- 1e5
lxj <- cumprod(c(l0,pj[-length(pj)]))
djk <- qjk * lxj

# step 3
lkxj <- rev(cumsum(rev(djk)))

# step 4
dj <- qj * lxj
Lj <- dj/Rj
nj <- diff(xj)
Ljk <- c(
  ((lkxj[-length(lkxj)] + lkxj[-1]) * nj)/2,
  (Djk[length(Djk)]/Dj[length(Dj)]) * Lj[length(Lj)]
  )

# step 5
Tkxj <- rev(cumsum(rev(Ljk)))

# step 6
ekxj <- Tkxj/lkxj

data.frame(group = xj, qjk = round(qjk,5), lkxj = round(lkxj), 
           djk = round(djk), Ljk = round(Ljk),
           Tkxj = round(Tkxj), ekxj = round(ekxj,2))

# function to create multiple decrement life table
mdlt <- function(group, all.deaths, cause.deaths, pop, radix=1e5){
  xj <- group
  
  # step 1
  Rj <- all.deaths/pop
  nj <- diff(group)
  qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
  Djk <- cause.deaths 
  Dj <- all.deaths
  qjk <- (Djk/Dj) * qj
  
  # step 2
  pj <- 1 - qj
  l0 <- radix
  lxj <- cumprod(c(l0,pj[-length(pj)]))
  djk <- qjk * lxj
  
  # step 3
  lkxj <- rev(cumsum(rev(djk)))
  
  # step 4
  dj <- qj * lxj
  Lj <- dj/Rj
  nj <- diff(xj)
  Ljk <- c(
    ((lkxj[-length(lkxj)] + lkxj[-1]) * nj)/2,
    (Djk[length(Djk)]/Dj[length(Dj)]) * Lj[length(Lj)]
  )
  
  # step 5
  Tkxj <- rev(cumsum(rev(Ljk)))
  
  # step 6
  ekxj <- Tkxj/lkxj
  
  data.frame(group = xj, qjk = round(qjk,5), lkxj = round(lkxj), 
             djk = round(djk), Ljk = round(Ljk),
             Tkxj = round(Tkxj), ekxj = round(ekxj,2))
}
# create multiple decrement life table for neoplasms
with(males, mdlt(group = age.group, all.deaths = deaths.all, 
                 cause.deaths = neoplasm.deaths, pop = pop))
mdlt.out <- with(males, 
                 mdlt(group = age.group, all.deaths = deaths.all,
                      cause.deaths = neoplasm.deaths, pop = pop))
# 27.17% of males born in Canada in 1991 predicted to die of neoplasm
mdlt.out[1,"Ljk"] / 1e5 * 100
# For an individual due to die of neoplasm, the life expectancy is at birth is
# 73.27 years
mdlt.out[1,"ekxj"]



# 13.3 --------------------------------------------------------------------

# cause-deleted life table


# function to cause-deleted life table

cdlt <- function(group, all.deaths, cause.deaths, pop, radix=1e5){
  Djk <- all.deaths - cause.deaths
  # annual death rate in population
  Rj <- (Djk/pop)
  # length of age group
  nj <- diff(group)
  # conditional probability of dying 
  qj <- c(nj*Rj[-length(Rj)] / (1 + (nj*Rj[-length(Rj)])/2), 1)
  # conditional probability of surviving
  pj <- 1 - qj
  # the radix: number of individuals in the OLT birth cohort
  l0 <- radix
  # number of survivors to age x
  lxj <- cumprod(c(l0,pj[-length(pj)]))
  # number of deaths in age group
  dj <- qj * lxj
  # number of person years in age group
  Lj <- dj/Rj
  # Number of person years after age x
  Txj <- sapply(1:length(Lj), function(x)sum(Lj[x:length(Lj)]))
  # Life expectancy at age x
  exj <- Txj/lxj
  # collect results into data frame and display
  data.frame(group = group, qj = round(qj,5), pj = round(pj, 5), 
             lxj = round(lxj), dj = round(dj), Lj = round(Lj), 
             Txj = round(Txj), exj = round(exj,2))
}


olt.out <- with(males, 
                olt(group = age.group, deaths = deaths.all, pop = pop))
circ.cdlt <- with(males, 
                  cdlt(group = age.group, all.deaths = deaths.all, 
                 cause.deaths = circ.deaths, pop = pop))
circ.mdlt <- with(males, 
                  mdlt(group = age.group, all.deaths = deaths.all, 
                 cause.deaths = circ.deaths, pop = pop))

# Table 13.8 summary indices of mortality - Circulatory
# Circulatory disease accounts for 40.09% of deaths
circ.mdlt[1,"lkxj"]/olt.out[1,"lxj"] * 100
# life expectancy (at birth) for those due to die of circulatory disease
circ.mdlt[1,"ekxj"]
# Eliminating circulatory diseases as a cause of death would increase overall
# life expectancy by 6.06 years
circ.cdlt[1,"exj"] - olt.out[1,"exj"]
# As a result of eliminating circulatory disease, the probability that a member
# of the MDLT cohort will survive to age 65 increases by about 13%
(circ.cdlt[circ.cdlt$group==65,"lxj"] - olt.out[olt.out$group==65,"lxj"]) / 
  circ.mdlt[1,"lkxj"] * 100
# Eliminating circulatory diseases as a cause of death would increase the life
# expectancy of those due to die of this cause by 15.12 years
(olt.out[1,"lxj"] * (circ.cdlt[1,"exj"] - olt.out[1,"exj"])) / circ.mdlt[1,"lkxj"]

