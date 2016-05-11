# Biostats in epidemiology
# ch 15 work


# 15.1 --------------------------------------------------------------------

# Logistic regression

# Example 15.1
str(breast.survival)
# This does not produce same results as Example 15.1
breast.survival$receptor.level <- relevel(breast.survival$receptor.level, ref = "high")
m1 <- glm(status ~ stage + receptor.level, data=breast.survival, family = binomial)
summary(m1)

# need to use table presented in Table 5.3
dat <- data.frame(t(as.data.frame(breast.stage)))
nms <- strsplit(row.names(dat),split = "\\.")
row.names(dat) <- NULL
dat$receptor.level <- factor(sapply(nms, function(x)x[1]), levels = c("high","low"))
dat$stage <- factor(sapply(nms, function(x)x[2]))

m1 <- glm(cbind(dead,alive)  ~ receptor.level + stage, data=dat, family = binomial)
summary(m1)
exp(confint(m1))
# estimated odds ratios
exp(coef(m1))
# Table 15.3
exp(coef(m1) + t(matrix(c(-1,1),nrow=2) %*% (qnorm(0.975) * summary(m1)$coefficients[,2])))

# prediction for receptor.level="low" and stage="II"
predict(m1, newdata = data.frame(receptor.level="low", stage="II"), type = "response")

# Model 15.7
m2 <- glm(cbind(dead,alive)  ~ receptor.level * stage, data=dat, family = binomial)
summary(m2)
# insignificant interaction implies conisderable evidence for homogeneity
mosaicplot(breast.stage)


# 15.2 Cox regression -----------------------------------------------------

library(survival)
str(breast.survival)
fit <- coxph(Surv(time, status) ~ receptor.level + stage, data = breast.survival) 
summary(fit)
