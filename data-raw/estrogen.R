# estrogen-endometrial cancer data
# Table 11.10, page 243

estrogen <- matrix(data = c(12,7,43,121),ncol=2,
                   dimnames = list("Case" = c("exposed","unexposed"),
                                   "Control" = c("exposed","unexposed")))

save(estrogen,file = "data/estrogen.rda")

estrogen2 <- matrix(data = c(1,0,10,1,10,1,10,1,2,0),ncol = 5,
                    dimnames = list("Case" = c("exposed","unexposed"),
                                    "Number of exposed controls" =0:4))

save(estrogen2,file = "data/estrogen2.rda")
