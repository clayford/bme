# oral contraceptives - myocardial infarction
# Table 11.2

oral <- matrix(c(29, 135, 205, 1607), ncol = 2, 
               dimnames = list("Myocardial infarction" = c("case","control"),
                               "Oral contraceptive" = c("yes","no")))
save(oral,file = "data/oral.rda")
