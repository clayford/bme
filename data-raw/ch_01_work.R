# probability function of X
X <- matrix(c(0,1,2,0.2,0.5,0.3),nrow = 3)

# expected value of X
crossprod(X[,1],X[,2])

# variance of X
varf <- function(x,m,p) (x-m)^2*p
sum(mapply(varf, x = X[,1], m = crossprod(X[,1],X[,2]), p = X[,2]))

# additive property of Normal dist

Zout.1 <- replicate(n = 1e5, sum(rnorm(10, mean = 1:10, sd = 1:10)))
Zout.2 <- replicate(n = 1e5, rnorm(1, mean = sum(1:10), sd = sqrt(sum((1:10)^2))))
op <- par(mfrow=c(2,1))
hist(Zout.1)
hist(Zout.2)
par(op)

# additive property of chi-square dist

Xout.1 <- replicate(n = 1e5, sum(rchisq(10, df = 1:10)))
Xout.2 <- replicate(n = 1e5, rchisq(1, df = sum(1:10)))
op <- par(mfrow=c(2,1))
hist(Xout.1)
hist(Xout.2)
par(op)

# Table 1.5
t1.5 <- cbind(dbinom(x = 0:10, size = 10, prob = 0.2),
      dbinom(x = 0:10, size = 20, prob = 0.1),
      dbinom(x = 0:10, size = 200, prob = 0.01),
      dpois(x = 0:10, lambda = 2))
apply(t1.5, 2, function(x)round(x*100,2))

# Figs 1.1(a) - 1.8(a)
plot(0:10, dbinom(x = 0:10, size = 10, prob = 0.3), type = "s", xlab = "Count", ylab = "Probability")
plot(0:10, dbinom(x = 0:10, size = 10, prob = 0.5), type = "s", xlab = "Count", ylab = "Probability")
plot(0:10, dbinom(x = 0:10, size = 100, prob = 0.03), type = "s", xlab = "Count", ylab = "Probability")
plot(0:15, dbinom(x = 0:15, size = 100, prob = 0.05), type = "s", xlab = "Count", ylab = "Probability")
plot(0:20, dbinom(x = 0:20, size = 100, prob = 0.1), type = "s", xlab = "Count", ylab = "Probability")

plot(0:10, dpois(x = 0:10, lambda = 3), type = "s", xlab = "Count", ylab = "Probability")
plot(0:15, dpois(x = 0:15, lambda = 5), type = "s", xlab = "Count", ylab = "Probability")
plot(0:20, dpois(x = 0:20, lambda = 10), type = "s", xlab = "Count", ylab = "Probability")

# Fig 1.9

p <- seq(0,1,0.01)
L <- p*(1 - p)^4
plot(p,L,type="l", ylim = c(0,0.1), xlab = expression(pi), ylab = expression(L(pi)), 
     las = 1)

