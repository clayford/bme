# Ovarian cancer data
# Table 9.6, p. 183
# data from a cohort study of women with stage II or stage IIIA ovarian cancer.
# grade is an indicator of the malignant potential of the tumor.
ovarian.cancer <- data.frame(grade = rep(c("Low","High"), c(15,20)), 
                      time = c(28, 89, 175, 195, 309, 377, 393, 421, 447, 462, 709, 744, 770,
                               1106, 1206,
                               34, 88, 137, 199, 280, 291, 299, 300, 309, 351, 358, 369, 369, 370,
                               375, 382, 392, 429, 451, 1119),
                      status = c(rep(1,5), rep(0,4), 1, rep(0,5),
                                 rep(1, 6), 0, 0, rep(1,9), 0, 1, 0))
# write.csv(ovarian, file = "ovarian.csv", row.names = FALSE)