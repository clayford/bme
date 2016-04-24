# breast cancer data

# Table 4.5(a), p. 99
breast <- matrix(c(23,25,31,113),ncol=2, dimnames = list(survival = c("dead","alive"), 
                                                         receptor.level=c("low","high")))
breast



# Table 5.3, p. 126
breast.stage <- array(data = c(2,10,5,50,
                               9,13,17,57,
                               12,2,9,6),
                      dim = c(2,2,3), dimnames = list(survival=c("dead","alive"),
                                                      receptor.level = c("low","high"),
                                                      stage=c("I","II","III")))

# Table 5.10, p. 140
# 
breast.receptor <- array(data = c(2,10,9,13,12,2,
                                  5,50,17,57,9,6),
                         dim = c(2,3,2), dimnames = list(survival=c("dead","alive"),
                                                         stage=c("I","II","III"),
                                                         receptor.level = c("low","high")))

# Table 9.1, p. 175
ILow <- c(50,51,51,53,53,54,54,55,56,56,57,60)
ILowc <- c(0,1,0,0,0,0,0,0,1,0,0,0)

IHigh <- c(10,34,34,47,47,49,49,rep(50,7),rep(51,6), rep(52,5), rep(53,6), rep(54,5),
           55,55,56,56,rep(57,5), rep(58,5), rep(59,4), 60, 60, 60)
IHighc <- c(1,1,0,1,1,rep(0, sum(c(2,7,6,5,6,5,2,2,5,5,4,3))))

IILow <- c(4,9,13,21,29,29,40,46,49,49,52,52,53,54,55,55,56,57,57,58,58,59,60)
IILowc <- c(0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
            1,0,0,0,0,0)

IIHigh <- c(11,16,21,23,23,24,33,33,36,36,36,37,45,46,49,49,rep(50,6),
            rep(51,4),rep(52,5), rep(53,5), rep(54,4), rep(55,4), rep(56,6), rep(57,4), 58, 58,
            rep(58,8), rep(59,5), rep(60,6))
IIHighc <- c(rep(1,10),0,1,1,1,rep(0, sum(c(2,6,4,5,5,4,4,6,4))), 1, 1,
             rep(0, sum(c(8,5,6))))

IIILow <- c(9, 12, 14, 15, 15, 17, 21, 22, 23, 23, 31, 34, 35, 53, 60)
IIILowc <- c(rep(1,4), 0, rep(1,8), 0, 0)

IIIHigh <- c(7, 9, 17, 21, 22, 22, 34, 34, 41, 49, 52, 55, 56, 58, 58, 59, 59)
IIIHighc <- c(0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, rep(0,5))

r0 <- c(12, 57, 23, 75, 15, 17)
breast.survival <- data.frame(time = c(ILow, IHigh, IILow, IIHigh, IIILow, IIIHigh),
                     status = c(ILowc, IHighc, IILowc, IIHighc, IIILowc, IIIHighc),
                     stage = c(rep("I", (12 + 57)), rep("II", (23 + 75)), rep("III", (15 + 17))),
                     receptor.level = rep(rep(c("low","high"), 3), r0))
breast.survival <- breast.survival[order(breast.survival$time),]
breast.survival$receptor.level <- relevel(breast.survival$receptor.level, ref = "low")
rm(list = ls(pattern = "I")); rm(r0)
# write.csv(tab9.1, file = "BreastCancer.csv", row.names = FALSE)