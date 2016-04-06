# Biostats in epidemiology
# ch 9 work

library(survival)

# Table 9.1
# ILow <- c(50,51,51,53,53,54,54,55,56,56,57,60)
# ILowc <- c(0,1,0,0,0,0,0,0,1,0,0,0)
# 
# IHigh <- c(10,34,34,47,47,49,49,rep(50,7),rep(51,6), rep(52,5), rep(53,6), rep(54,5),
#            55,55,56,56,rep(57,5), rep(58,5), rep(59,4), 60, 60, 60)
# IHighc <- c(1,1,0,1,1,rep(0, sum(c(2,7,6,5,6,5,2,2,5,5,4,3))))
# 
# IILow <- c(4,9,13,21,29,29,40,46,49,49,52,52,53,54,55,55,56,57,57,58,58,59,60)
# IILowc <- c(0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,
#             1,0,0,0,0,0)
# 
# IIHigh <- c(11,16,21,23,23,24,33,33,36,36,36,37,45,46,49,49,rep(50,6),
#             rep(51,4),rep(52,5), rep(53,5), rep(54,4), rep(55,4), rep(56,6), rep(57,4), 58, 58,
#             rep(58,8), rep(59,5), rep(60,6))
# IIHighc <- c(rep(1,10),0,1,1,1,rep(0, sum(c(2,6,4,5,5,4,4,6,4))), 1, 1, 
#              rep(0, sum(c(8,5,6))))
# 
# IIILow <- c(9, 12, 14, 15, 15, 17, 21, 22, 23, 23, 31, 34, 35, 53, 60)
# IIILowc <- c(rep(1,4), 0, rep(1,8), 0, 0)
# 
# IIIHigh <- c(7, 9, 17, 21, 22, 22, 34, 34, 41, 49, 52, 55, 56, 58, 58, 59, 59)
# IIIHighc <- c(0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, rep(0,5))
# 
# r0 <- c(12, 57, 23, 75, 15, 17)
# tab9.1 <- data.frame(time = c(ILow, IHigh, IILow, IIHigh, IIILow, IIIHigh),
#                      status = c(ILowc, IHighc, IILowc, IIHighc, IIILowc, IIIHighc),
#                      stage = c(rep("I", (12 + 57)), rep("II", (23 + 75)), rep("III", (15 + 17))),
#                      receptor = rep(rep(c("Low","High"), 3), r0))
# 
# write.csv(tab9.1, file = "BreastCancer.csv", row.names = FALSE)
tab9.1 <- read.csv("BreastCancer.csv")

# Table 9.2
# create a Surv object
SurvBC <- Surv(tab9.1$time, tab9.1$status)
surv.all <- survfit(SurvBC ~ 1, conf.type = "log-log")
summary(surv.all)

# Figure 9.2
plot(surv.all, ylim=c(0.5,1.0), ylab = "Survival probability", xlab = "Time")
