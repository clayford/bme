# oral contraceptives - myocardial infarction
# Table 11.2

oral <- matrix(c(29, 135, 205, 1607), ncol = 2, 
               dimnames = list("Myocardial infarction" = c("case","control"),
                               "Oral contraceptive" = c("yes","no")))
save(oral,file = "data/oral.rda")

# Table 11.5(a)

oral.age <- array(data = c(13,95,14,614,10,35,98,692,6,5,93,301), dim = c(2,2,3), 
                  dimnames = list("disease" = c("case","control"),
                                  "OC" = c("yes","no"),
                                  "age group" = c("25-34","35-44","45-49")))
save(oral.age,file = "data/oral.age.rda")
