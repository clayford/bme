# ch 2
library(epitools)

# 2.3.3
# Table 2.2(a) - 2.2(e)

tab.a <- expand.grid(D = 1:2, E = 1:2, F = 1:2)
# fill by column, going across tables, left to right
tab.a$Freq <- c(70, 30, 40, 60, 140, 60, 80, 120)
t.out <- xtabs(Freq ~ D + E + F, data = tab.a)
xtabs(Freq ~ D + E, data = tab.a) # F = . (crude table)

# effect measures
riskratio(t.out[,,1])$measure[2,1]
riskratio(t.out[,,2])$measure[2,1]
oddsratio(t.out[,,1])$measure[2,1]
oddsratio(t.out[,,2])$measure[2,1]

# risk difference
p.out <- prop.table(t.out[,,1], margin = 2)
p.out[1,1] - p.out[1,2]

RD <- function(t){
  p.out <- prop.table(t, margin = 2)
  p.out[1,1] - p.out[1,2]
}

RD(t = t.out[,,1])

effMeas <- function(tab){
  # need to transpose tables and reverse cols and rows to use with epitools functions
  t1 <- tab[,,1]
  t2 <- tab[,,2]
  ct <- margin.table(tab, margin = c(1,2))
  list("F=1" = tab[,,1],"RD" = RD(t1),
                        "RR" = riskratio(t(t1), rev = "b")$measure[2,1],
                        "OR" = oddsratio(t(t1), rev = "b")$measure[2,1],
        "F=2" = tab[,,2], "RD" = RD(t2),
                         "RR" = riskratio(t(t2), rev = "b")$measure[2,1],
                         "OR" = oddsratio(t(t2), rev = "b")$measure[2,1],
        "F=." = ct, "RD" = RD(ct),
                    "RR" = riskratio(t(ct), rev = "b")$measure[2,1],
                    "OR" = oddsratio(t(ct), rev = "b")$measure[2,1])
}

# Generate tables 2.2(a) - 2.2(e)

tab <- expand.grid(D = 1:2, E = 1:2, F = 1:2)
# fill by column, going across tables, left to right

# a
tab$Freq <- c(70, 30, 40, 60, 140, 60, 80, 120)
t.out <- xtabs(Freq ~ D + E + F, data = tab)
effMeas(t.out)
# b
tab$Freq <- c(70, 30, 40, 60, 160, 40, 80, 120)
t.out <- xtabs(Freq ~ D + E + F, data = tab)
effMeas(t.out)
# c
tab$Freq <- c(70, 30, 80, 120, 160, 40, 40, 60)
t.out <- xtabs(Freq ~ D + E + F, data = tab)
effMeas(t.out)
# d
tab$Freq <- c(90, 10, 60, 40, 80, 120, 20, 180)
t.out <- xtabs(Freq ~ D + E + F, data = tab)
effMeas(t.out)
# e
tab$Freq <- c(90, 10, 120, 80, 30, 170, 10, 90)
t.out <- xtabs(Freq ~ D + E + F, data = tab)
effMeas(t.out)
