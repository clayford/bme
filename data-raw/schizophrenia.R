# Table 12.2(a)

age.group <- c(paste(seq(10,70,10),seq(10,70,10) + 9,sep = "-"),"80+")
cohort.deaths <- c(2,55,32,21,27,19,25,9)
cohort.py <- c(285.1, 4179.1, 3291.2, 1994.7, 1498.9, 763.5, 254.4, 46.7)
alberta.deaths <- c(267, 421, 306, 431, 836, 1364, 1861, 1797)
alberta.pop <- c(201825, 263175, 176140, 114715, 93315, 60835, 34250, 12990)
schizophrenia <- data.frame(age.group, cohort.deaths, cohort.py, alberta.deaths, alberta.pop)

save(schizophrenia, file="data/schizophrenia.rda")
