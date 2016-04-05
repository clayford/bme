# Biostats in epidemi0logy
# ch 3 work

# Table 3.1
a <- 0:10
data.frame(
  a,
  "p(A = a)" = round(dbinom(x = a, size = 10, prob = 0.4),4)*100,
  "p(A <= a)" = round(pbinom(q = a, size = 10, prob = 0.4),4)*100,
  "p(A >= a)" = round(
      dbinom(x = a, size = 10, prob = 0.4) +
      pbinom(q = a, size = 10, prob = 0.4, lower.tail = FALSE),
      4)*100,
  check.names = FALSE
)

# Example 3.2
# Test null that p = 0.4
binom.test(x = 2, n = 10, p = 0.4)

# Example 3.3
prop.test(x = 10, n = 50, correct = FALSE)
# "The confidence interval is computed by inverting the score test."
# This appears to be the same as the "implicit method" in the book

# The Wald method (or explicit method) from the book:
moe <- qnorm(0.975)*sqrt((0.2*0.8)/50)
0.2 + c(-moe,moe)

# Example 3.4
p.out <- prop.test(x = 2, n = 10, correct = FALSE)
round(p.out$conf.int*100,2)

moe <- qnorm(0.975)*sqrt((0.2*0.8)/10)
0.2 + c(-moe,moe)

# Example 3.5
binom.test(x = 2, n = 10, p = 0.4)$p.value
prop.test(x = 2, n = 10, p = 0.4, correct = FALSE)$p.value

mapply(binom.test, x=c(2,5,10), n=c(10,25,50), p=0.4)["p.value",]
mapply(prop.test, x=c(2,5,10), n=c(10,25,50), p=0.4, correct=FALSE)["p.value",]


mapply(prop.test, x=c(2,5,10), n=c(10,25,50), p=0.4, correct=FALSE)

m.out <- mapply(prop.test, x=c(2,5,10), n=c(10,25,50), p=0.4, correct=FALSE)
str(m.out)

pv <- function(x, n){
  prop.test(x, n, p=0.4, correct=FALSE)$p.value
}
mapply(pv, x=c(2,5,10), n=c(10,25,50))


# Example 3.6
prop.test(x = 54, n = 192, correct = FALSE)

library(epitools)


r243 <- matrix(c(12,2,7,9), 2, 2)
dimnames(r243) <- list(Diarrhea = c("Yes", "No"),
                       "Antibody level" = c("Low", "High")
)
r243
r243b <- t(r243)
r243b
epitab(r243, rev = "b", verbose = TRUE)
epitab(r243, method="riskratio",rev = "b", verbose = TRUE)
epitab(matrix(c(41, 15, 28010, 19017),2,2)[2:1,],
       method="rateratio", verbose = TRUE)