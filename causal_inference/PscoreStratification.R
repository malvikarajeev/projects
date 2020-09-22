library(car)
library(ATE)
## function for SRE
## estimation
Neyman_SRE = function(z, y, x)
{
  xlevels = unique(x)
  K       = length(xlevels)
  PiK     = rep(0, K)
  TauK    = rep(0, K)
  varK    = rep(0, K)
  for(k in 1:K)
  {
    xk         = xlevels[k]
    zk         = z[x == xk]
    yk         = y[x == xk]
    PiK[k]     = length(zk)/length(z)
    TauK[k]    = mean(yk[zk==1]) - mean(yk[zk==0])
    varK[k]    = var(yk[zk==1])/sum(zk) + 
      var(yk[zk==0])/sum(1 - zk)
  }
  
  return(c(sum(PiK*TauK), sqrt(sum(PiK^2*varK))))
}


data(nhanes_bmi)
z = nhanes_bmi$School_meal
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])
x = scale(x)

## simple regression analysis
diffinmeans = lm(y ~ z)
fisher.ancova = lm(y ~ z + x)
lin.ancova = lm(y ~ z + x + z*x)

res.regression = c(coef(diffinmeans)[2], hccm(diffinmeans)[2, 2]^0.5,
                   coef(fisher.ancova)[2], hccm(fisher.ancova)[2, 2]^0.5,
                   coef(lin.ancova)[2], hccm(lin.ancova)[2, 2]^0.5)
res.regression = matrix(res.regression,
                        nrow = 2, ncol = 3)
rownames(res.regression) = c("est", "se")
colnames(res.regression) = c("naive", "fisher", "lin")
round(res.regression, 2)                  

## propensity score stratification
pscore = glm(z ~ x, family = binomial)$fitted.values
par(mfrow = c(1, 3))
hist(pscore[z==1], freq = FALSE, col = "grey",
     border = NA, xlab = "estimated propensity score",
     ylab = "", yaxt = "n", breaks = 30, 
     main = "breaks = 30", 
     xlim = c(0, 1), ylim = c(0, 5))
hist(pscore[z==0], freq = FALSE, 
     add = TRUE, breaks = 30)

## five strata with cutoff points (1:5)/5
hist(pscore[z==1], freq = FALSE, col = "grey",
     border = NA, xlab = "estimated propensity score",
     ylab = "", yaxt = "n", breaks = 5, 
     main = "breaks = 5", xlim = c(0, 1))
hist(pscore[z==0], freq = FALSE, 
     add = TRUE, breaks = 5)

## ten strata with cutoff points (1:10)/10
hist(pscore[z==1], freq = FALSE, col = "grey",
     border = NA, xlab = "estimated propensity score",
     ylab = "", yaxt = "n", breaks = 10, 
     main = "breaks = 10", 
     xlim = c(0, 1), ylim = c(0, 3))
hist(pscore[z==0], freq = FALSE, 
     add = TRUE, breaks = 10)

pscore.strata = cut(pscore, breaks = (0:5)/5, 
                    labels = 1:5)
Neyman_SRE(z, y, pscore.strata)

pscore.strata = cut(pscore, breaks = (1:10)/10, 
                    labels = 1:9)
Neyman_SRE(z, y, pscore.strata)

## discretized by the quantiles of the pscore
n.strata = c(5, 10, 20, 50, 80)
strat.res = sapply(n.strata, 
                   FUN = function(nn){
                     q.pscore = quantile(pscore, (1:(nn-1))/nn)
                     ps.strata = cut(pscore, breaks = c(0,q.pscore,1), 
                                     labels = 1:nn)
                     Neyman_SRE(z, y, ps.strata)             
                   })
rownames(strat.res) = c("est", "se")
colnames(strat.res) = n.strata
round(strat.res, 2)

## balance check for K=5
Bcheck = sapply(1:dim(x)[2],
                FUN = function(px){
                  q.pscore = quantile(pscore, (1:4)/5)
                  ps.strata = cut(pscore, breaks = c(0,q.pscore,1), 
                                  labels = 1:5)
                  Neyman_SRE(z, x[, px], ps.strata)
                })
rownames(Bcheck) = c("est", "se")
colnames(Bcheck) = 1:dim(x)[2]
round(Bcheck, 2)

                     
