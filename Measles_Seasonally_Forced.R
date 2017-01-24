
#Required Packages
library(deSolve)
require(pracma)
require(FME)
require(stats4)
require(bbmle)

#Data reading
Data = read.table(file.choose())
summary(Data)
head(Data)
Data = Data[1:469,2]/5.1e+07

#Reporting Function
report=function(out, gamma = gamma){
  I= out[,4]
  C=c()
  for (i in 0:(468+1456) ){
    C[i+1]=0
    for (j in 1:7) {
      x=(7*i)+j
      C[i+1]=C[i+1] + gamma*I[x]*(1/(365.25*7))  
    }
  }
  return(C)
}


"MODEL"
pars <- c(mub=1/60, beta= 600, mud = 1/60, lambda = 365/8, gamma = 365/14,
          va = 0, s = 1/23, e = 0, i = 1e-4, r = 1-1/23-1e-4 , beta2 = .2, phi = .5)


pars <- as.vector(pars)



  pars = as.vector(pars)
  
  #Initial Conditions
  s = as.vector(pars[7])
  e = as.vector(pars[8])
  i = as.vector(pars[9])
  r = as.vector(pars[10])
  n = s+e+i+r 
  
  state <- c(S = s , E = e , I = i, R = r, N = n)
  
  
  sir_rhs=function(t,state,pars){
    
    
    with(as.list(c(state, pars)),{
      #rates of change
      dS <-  pars[1]*(1-pars[6]) -  (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12]))))  - pars[3]*S
      dE <- (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12])))) - (pars[4]+pars[3])*E
      dI <-  pars[4]*E - (pars[5]+pars[3])*I
      dR <- pars[5]*I + (pars[1]*pars[6]) - pars[3]*R
      dN <- pars[1] - pars[3]*N  
      
      # return the rate of change
      return(list(c(dS, dE, dI, dR, dN)))
    })
  }
o = 37/((469+1456)*7)          
times <- seq(0, 37, by = o)
times <- times[1:((469+1456)*7)] # 37 years in days  (length = 13475)

#Realisation
outt<- ode(y = state, times = times, func = sir_rhs, parms = pars)

"Cost Function"
least_squares=function(parameter, data){
  pars <- as.vector(parameter)
  
  mub=pars[1]
  beta=pars[2]
  mu = pars[3] 
  lambda = pars[4]
  gamma = pars[5] 
  va =pars[6] 
  beta2 = pars[11]
  phi = pars[12]
  s = pars[7]
  e = pars[8]
  i = pars[9]
  r = pars[10]
  n = s+e+i+r
  
  
  state <- c(S = s , E = e , I = i, R = r, N = n)
  
  #Model representation
  
  sir_rhs=function(t,state,pars){
    
    
    with(as.list(c(state, pars)),{
      #rates of change
      dS <-  pars[1]*(1-pars[6]) -  (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12]))))  - pars[3]*S
      dE <- (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12])))) - (pars[4]+pars[3])*E
      dI <-  pars[4]*E - (pars[5]+pars[3])*I
      dR <- pars[5]*I + (pars[1]*pars[6]) - pars[3]*R
      dN <- pars[1] - pars[3]*N  
      
      # return the rate of change
      return(list(c(dS, dE, dI, dR, dN)))
    })
  }
  
  out <- ode(y = state, times = times, func = sir_rhs, parms = pars)
  
  OUT<- report(out,gamma)
  OUT = OUT[-(1:1456)]
  data = Data
  ESS = 0
  error = data - OUT #assuming out[,2] is for incidence
  ESS = sum(error^2)
  return(ESS)
}

least_squares(pars,Data)

# I initialy use the next set of wriiten codes for my least sqaure minimization, However
#Using the fminsearch() of the pracma package doesnt allow putting putting a bound on the parameters 
# Does under the Sensitivity section is a much better parameter estimation
"Test Minimizers Done"


require(pracma)
Model_fit=function(pars){
  estpars<- fminsearch(least_squares, pars, minimize = TRUE, dfree = F,
                       maxiter = 10000, tol = .Machine$double.eps^(2/3))
  return(estpars)
}

Newpars = Model_fit(pars)

Error=function(X,Y){
  ESS = 0
  error = X - Y #assuming out[,2] is for incidence
  ESS = sum(error^2)
  return(ESS)
}

parest = function (pars){
  pars = pars
  estpars<-Model_fit(pars)$xval
  n = 1
  while ( (Error(least_squares(pars),least_squares(estpars))   > 0.05 * least_squares(pars))|(n<=10)){
    pars<-estpars
    estpars<-Model_fit(pars)$xval
    n = n + 1
  }
  return(estpars)
}

newp<-parest(pars)


############################################################################

###########       Maximum Likelihood Estimation       #############

############################################################################



log_likelihood=function(mub,beta,mud,lambda,gamma,va,s,e,i,r, beta2, phi, mu, sigma){
  
  pars <- c()
 
  pars[5] = gamma
  pars[1] = mub
  pars[2] = beta
  pars[3] = mud
  pars[4] = lambda
  gamma = pars[5] 
  pars[6] = va
  pars[11] = beta2 
  pars[12] = phi 
  pars[7] =  s 
  pars[8] = e
  pars[9] = i
  pars[10] = r
  n = s+e+i+r
  
  
  state <- c(S = s , E = e , I = i, R = r, N = n)
  
  #Model representation
  
  sir_rhs=function(t,state,pars){
    
    
    with(as.list(c(state, pars)),{
      #rates of change
      dS <-  pars[1]*(1-pars[6]) -  (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12]))))  - pars[3]*S
      dE <- (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12])))) - (pars[4]+pars[3])*E
      dI <-  pars[4]*E - (pars[5]+pars[3])*I
      dR <- pars[5]*I + (pars[1]*pars[6]) - pars[3]*R
      dN <- pars[1] - pars[3]*N  
      
      # return the rate of change
      return(list(c(dS, dE, dI, dR, dN)))
    })
  }
  
  times <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]
  
  out <- ode(y = state, times = times, func = sir_rhs, parms = pars)
  
  OUT<- report(out,gamma)
  OUT = OUT[-(1:1456)]
  data = Data
  Resd =  data - OUT #Residuals assumed with mean 0
  R = suppressWarnings(dnorm(Resd , mu, sigma, log = T))
  -sum(R)
}

fit = mle(log_likelihood, start = list(mub=1/60, beta= 600, mud = 1/60, lambda = 365/8, gamma = 365/14,
                                       va = 0, s = 1/23, e = 0, i = 1e-4, r = 1-1/23-1e-4 , beta2 = .2, phi = .5,  mu = 0, sigma = 3.637182)
          , fixed = list(gamma = 365/14, mu = 0, va = 0), lower = rep(0,12),method = "L-BFGS-B" , nobs = 469)

fit
summary(fit)
AIC(fit)
BIC(fit)
logLik(fit)

















"Sensitivity Analysis"
########################

gamma = 365/14

pars <- c(mub=1/60, beta= 600, mud = 1/60, lambda = 365/8, gamma = 365/14,
          va = 0, s = 1/23, e = 0, i = 1e-4, r = 1-1/23-1e-4 , beta2 = .2, phi = .5)
#SEIRR Gives the Model realisation weekly and appends the Dataframe with Observation (L)
SEIRR<-function(pars){
  pars = as.vector(pars)
  #Initial Conditions
  s = as.vector(pars[7])
  e = as.vector(pars[8])
  i = as.vector(pars[9])
  r = as.vector(pars[10])
  n = s+e+i+r 
  
  state <- c(S = s , E = e , I = i, R = r, N = n)
  
  sir_rhs=function(t,state,pars){
    
    
    with(as.list(c(state, pars)),{
      #rates of change
      dS <-  pars[1]*(1-pars[6]) -  (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12]))))  - pars[3]*S
      dE <- (S*I*pars[2]*(1 + pars[11]*cos(2*pi*(t - pars[12])))) - (pars[4]+pars[3])*E
      dI <-  pars[4]*E - (pars[5]+pars[3])*I
      dR <- pars[5]*I + (pars[1]*pars[6]) - pars[3]*R
      dN <- pars[1] - pars[3]*N  
      
      # return the rate of change
      return(list(c(dS, dE, dI, dR, dN)))
    })
  }

  
  times <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]

  out <- ode(y = state, times = times, func = sir_rhs, parms = pars)
  
  new=as.data.frame(out)
  
  
  time <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]
  timesD<- cbind(t=0:13474,time) 
  ii <- which (timesD[,1] %in% seq(0, 13474, by = 7))
  
  Dat <- cbind(new[ii,], L=report(out,gamma))[-(1:1456),]
  return(Dat)
}

out <- SEIRR(pars)

#time Dataset


time <- seq(0, 37, by = 37/((469+1456)*7) )[1:((469+1456)*7)]
timesD<- cbind(t=0:13474,time) 
ii <- which (timesD[,1] %in% seq(0, 13474, by = 7))[-(1:1456)]


DataT <- cbind(time = timesD[ii,2], L= Data) # Appending time and Data


# Model Cost
SEIRRcost<- function(pars){
  out<-SEIRR(pars)
  cost <- modCost(model = out, obs = DataT)
  return(cost)
}


#Parameter Estimation
Fit <- modFit(f = SEIRRcost, p = pars, lower = rep(0,12), upper = c(rep(Inf,10),1,1)  )


# Local Sensitivities
Sfun<-sensFun(func = SEIRR, parms = pars, senspar = c(1:6,11,12))
summary(Sfun)
plot(Sfun, which = c("L","S","E","I","R","N"), xlab ="time", lwd = 2)

# Identifiability of parameters
ident <- collin(Sfun)
ident
plot(ident, log = "y")

## 20 = magical number above which there are identifiability problems
abline(h = 20, col = "red")

#Manual Iteration for parameter estimates
pars = Fit$par
gamma = pars[5]
Fit <- modFit(f = SEIRRcost, p = pars, lower = rep(0,12), upper = c(rep(Inf,10),1,1) )


# Methods can be any of these c("Marq", "Port", "Newton","Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Pseudo", "bobyqa")