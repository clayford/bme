# Biostats in epidemiology
# ch 6 work

#Risk ratio methods for closed cohort data
library(epitools)
t45a <- matrix(c(23,25,31,113),ncol=2, dimnames = list(survival = c("dead","alive"), 
                                                      receptor.level=c("low","high")))
# example 6.1
# the following estimates risk ratio and provides CI
riskratio(t(t45a), rev = "both")

# risk ratio by hand
a <- breast[1,]
r <- apply(breast,2,sum)
m <- apply(breast,1,sum)
RR <- (a[1]*r[2])/(a[2]*r[1])


# Wald test of association using epitools
# RR <- riskratio(t(t45a), rev = "both")$measure[2,1]
# r1 <- apply(t45a,2,sum)[1]
# r2 <- apply(t45a,2,sum)[2]
# m1 <- apply(t45a,1,sum)[1]
# m2 <- apply(t45a,1,sum)[2]

X_w <- ((log(RR)^2) * r[1] * r[2] * m[1]) / (sum(t45a) * m[2])
pchisq(X_w, df = 1, lower.tail = FALSE)

# LRT of association
X_lr <- 2*sum(t45a*log(t45a/epitools::expected(t45a)))
pchisq(X_lr, df = 1, lower.tail = FALSE)

# function for this?
risk.ratio.test <- function(data, conf.level=0.95, Wald=TRUE){
  dname <- deparse(substitute(data))
  alternative <- "two.sided"

    a <- data[1,]
    b <- data[2,]
    r <- apply(data,2,sum)
    m <- apply(data,1,sum)
    # risk ratio
    if(a[1]==0 || a[2]==0){
      est <- ((a[1] + 0.5) *r[2])/((a[2] + 0.5)*r[1])
    } else {
      est <- (a[1]*r[2])/(a[2]*r[1])  
    }
    names(est) <- "risk ratio"
    null <- 1
    names(null) <- names(est)
    # variance
    v <- b[1]/(a[1]*r[1]) + b[2]/(a[2]*r[2])
    alpha <- (1-conf.level)/2
    CINT <- exp(log(est) + c(-1,1) * qnorm(1 - alpha) * sqrt(v))
    attr(CINT, "conf.level") <- conf.level
    if(Wald){
      STATISTIC <- ((log(est)^2) * r[1] * r[2] * m[1]) / (sum(data) * m[2])
      p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  } else {
    # LRT of association
    STATISTIC <- 2*sum(data * log(data/epitools::expected(data)))
    p.value <- pchisq(STATISTIC, df = 1, lower.tail = FALSE)
  }
  names(STATISTIC) <- "X-squared"
  METHOD <- paste(if(Wald) "Wald" else "Likelihood Ratio", "Test of association for risk ratio")
  RVAL <- list(statistic = STATISTIC, parameter = c(df = 1), p.value = p.value, estimate = est, 
               null.value = null,
               conf.int = CINT, alternative = alternative,
               method = METHOD, 
               data.name = dname) 
  class(RVAL) <- "htest"
  return(RVAL)
}
# Example 6.1
risk.ratio.test(breast)
risk.ratio.test(breast, Wald = FALSE)


# 6.3 ---------------------------------------------------------------------

# Mantel-Haenszel estimate of the risk ratio

# table 5.3
# receptor level breast cancer
t53 <- array(data = c(2,10,5,50,
                      9,13,17,57,
                      12,2,9,6),
             dim = c(2,2,3), dimnames = list(survival=c("dead","alive"),
                                             ReceptorLevel = c("low","high"),
                                             stage=c("I","II","III")))
# Mantel-Haenszel estimate of the risk ratio
R <- sum(t53[1,1,] * apply(t53[,2,],2,sum) / apply(t53,3,sum))
S <- sum(t53[1,2,] * apply(t53[,1,],2,sum) / apply(t53,3,sum))
RR_mh <- R/S

# variance of log(RR)
Tjsum <- sum(
  ((apply(t53[,2,],2,sum) * apply(t53[,1,],2,sum) * apply(t53[1,,],2,sum)) - 
  (t53[1,1,] * t53[1,2,] * apply(t53,3,sum))
) / apply(t53,3,sum)^2
)

vRR_mh <- Tjsum / (R * S)

# 95% CI
exp(log(RR_mh) + c(-1,1) * qnorm(0.975) * sqrt(vRR_mh))

# function for MH est of RR with 95% CI
# disease on rows (yes in first row); exposure on columns; stratum is 3rd layer 
mh.rr <- function(x){
  if(length(dim(x)) != 3 || dim(x)[-3] != c(2,2)) stop("This function is for 2 x 2 x K tables")
  # Mantel-Haenszel estimate of the risk ratio
  R <- sum(x[1,1,] * apply(x[,2,],2,sum) / apply(x,3,sum))
  S <- sum(x[1,2,] * apply(x[,1,],2,sum) / apply(x,3,sum))
  RR_mh <- R/S
  
  # variance of log(RR)
  Tjsum <- sum(
    ((apply(x[,2,],2,sum) * apply(x[,1,],2,sum) * apply(x[1,,],2,sum)) - 
       (x[1,1,] * x[1,2,] * apply(x,3,sum))
    ) / apply(x,3,sum)^2
  )
  
  vRR_mh <- Tjsum / (R * S)
  list(MH.risk.ratio = RR_mh, 
       CI = exp(log(RR_mh) + c(-1,1) * qnorm(0.975) * sqrt(vRR_mh)))
}
mh.rr(t53)

# Using the Hmisc function mhgr
# need to transform data to one record per obs
t53df <- as.data.frame(as.table(t53))
surv <- rep(rep(c(1,0),6), t53df$Freq)
t53df2 <- t53df[rep(1:12, t53df$Freq),]
t53df2$surv <- surv
t53df2$Freq <- NULL

with(t53df2, Hmisc::mhgr(y = surv, group = ReceptorLevel, strata = stage))

