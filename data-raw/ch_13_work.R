# Biostats in epidemiology
# ch 13 work

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

