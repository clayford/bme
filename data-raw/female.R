# female death rates
# Table 12.3
# page 258


females <- matrix(c(1.4, 0.9, 0.9, 0.7, 0.6,
                    3.2, 2.1, 2.1, 1.6, 1.2,
                    6.6, 5.3, 4.9, 4.2, 3.4,
                    16.1, 13.4, 11.2, 9.8, 8.4,
                    42.8, 35.1, 29.4, 24.9, 21.5),
                  ncol=5)
dimnames(females) <- list("Year" = seq(1950, 1990, 10),
                          "Age Group" = c("30-34","40-44", "50-54", "60-64","70-74"))
save(females, file = "data/females.rda")
