#Maximum likelihood fitting
x = 1:100
y = 3*x + 5 + rnorm(100, 0, 4)

lmod = lm(y~x)

summary(lmod)

#Data fitting

ll = function(b0, b1, mu, sigma){
  Resd = y - (b0 + b1*x)
  R = suppressWarnings(dnorm(Resd , mu, sigma, log = T))
  
  -sum(R)
}

require(stats4)
fit = fit = mle(ll, start = list(b0 = 5.135368, b1 = 3.637182 , mu = 0, sigma = 3.637182), fixed = list(mu = 0), nobs = length(y))

fit
summary(fit)
AIC(fit)
BIC(fit)
logLik(fit)
