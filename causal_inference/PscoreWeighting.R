ipw.est = function(z, y, x, truncpscore = c(0, 1))
{
  ## fitted propensity score
  pscore   = glm(z ~ x, family = binomial)$fitted.values
  pscore   = pmax(truncpscore[1], pmin(truncpscore[2], pscore))
  
  ace.ipw0 = mean(z*y/pscore - (1 - z)*y/(1 - pscore))
  ace.ipw  = mean(z*y/pscore)/mean(z/pscore) - 
    mean((1 - z)*y/(1 - pscore))/mean((1 - z)/(1 - pscore))

  return(c(ace.ipw0, ace.ipw))     
}


ipw.boot = function(z, y, x, n.boot = 200, truncpscore = c(0, 1))
{
  point.est  = ipw.est(z, y, x, truncpscore)
  
  ## nonparametric bootstrap
  n.sample   = length(z)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot, 
                         {id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                         ipw.est(z[id.boot], y[id.boot], x[id.boot, ])})
  boot.var   = apply(boot.est, 1, sd)
  
  res = rbind(point.est, boot.var)
  rownames(res) = c("est", "se")
  colnames(res) = c("HT", "Hajek")
  
  return(res)
}


library(ATE)
data(nhanes_bmi)
z = nhanes_bmi$School_meal
y = nhanes_bmi$BMI
x = as.matrix(nhanes_bmi[, -c(1, 2)])
x = scale(x)

ipw.boot(z, y, x)
ipw.boot(z, y, x, truncpscore = c(0.01, 0.99))
ipw.boot(z, y, x, truncpscore = c(0.05, 0.95))
ipw.boot(z, y, x, truncpscore = c(0.1, 0.9))

## balance check based on Hajek
Bcheck = sapply(1:dim(x)[2],
                FUN = function(px){
                  ipw.boot(z, x[, px], x)[, 2]
                })
rownames(Bcheck) = c("est", "se")
colnames(Bcheck) = 1:dim(x)[2]
round(Bcheck, 2)

