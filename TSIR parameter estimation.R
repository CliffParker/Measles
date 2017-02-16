#Packages
require(pracma)
#Cost function
parcost = function(alpha1 = 1, zeta =1.573e-6, Z0 = -11259.898, I0 =11389.669, rmod = rep(1.349859,26),dta = Data){

  #Creating R vector
  r = rep(NA,226)
  r[1:226] = rep( rmod , 9)[1:226]
  #Setting initial values
  dta[1,10] = Z0
  dta[1,11] = I0
  
  Y = log(dta$I)
  Ycap = log(r) + alpha1 * log(dta$`It-1`) + zeta * dta$`Zt-1`  # without error
  
  ESS = 0
  error = Y - Ycap 
  ESS = sum(error^2)
  return(ESS)
  
}

#Parameter estimation section
Model_fit=function(par = c( alpha1 = 1, 
                            zeta = 1.573e-6, Z0 = -11259.898, I0 = 11389.669, rmod = rep(1.349859,26))){
  estpars<- fminsearch(parcost, pars <- par, 
            minimize = T, dfree = F,
                       maxiter = 1e+6, tol = .Machine$double.xmin)
  return(estpars)
}

#Iterative parameter estimation
Fitit = function(n){
  est = Model_fit()$xval
  for(i in 1:n){
    est = Model_fit(est)$xval
  }
  return(Model_fit(est))
}




log_lik = function(alpha1 = 1, zeta =1.573e-6, Z0 = -11259.898, I0 =11389.669,
                   rmod1,rmod2,rmod3,rmod4,rmod5,rmod6,rmod7,rmod8,rmod9,rmod10,rmod11,rmod12,rmod13,rmod14,rmod15,
                   rmod16,rmod17,rmod18,rmod19,rmod20,rmod21,rmod22,rmod23,rmod24,rmod25,rmod26, mu = 0, sigma = 57499.67){
  
  r = rep(NA,226)
  rmod = c(rmod1,rmod2,rmod3,rmod4,rmod5,rmod6,rmod7,rmod8,rmod9,rmod10,rmod11,rmod12,rmod13,rmod14,rmod15,
           rmod16,rmod17,rmod18,rmod19,rmod20,rmod21,rmod22,rmod23,rmod24,rmod25,rmod26 )
  r[1:226] = rep( rmod , 9)[1:226]
  dta =Data
  dta[1,10] = Z0
  dta[1,11] = I0
  
  Y = log(dta$I)
  Ycap = log(r) + alpha1 * log(dta$`It-1`) + zeta * dta$`Zt-1`  # without error
  
  Resd =  Y - Ycap #Residuals assumed with mean 0
  R = suppressWarnings(dnorm(Resd , mu, sigma, log = T))
  -sum(R)
  
}
require(bbmle)
fit = mle2(log_lik, start = list(alpha1 = 1, zeta =1.573e-6, Z0 = -11259.898, I0 =11389.669, 
                                 rmod1 = 1.2,rmod2= 1.3,rmod3= 1.33,rmod4= 1.35,rmod5= 1.3,rmod6= 1.3,rmod7= 1.3,rmod8= 1.3
                                 ,rmod9= 1.3,rmod10= 1.3,rmod11= 1.3,rmod12= 1.3,rmod13= 1.3,rmod14= 1.3,rmod15= 1.3,
                                 rmod16= 1.3,rmod17= 1.3,rmod18= 1.3,rmod19= 1.3,rmod20= 1.3,rmod21= 1.3,rmod22= 1.3
                                 ,rmod23= 1.3,rmod24= 1.3,rmod25= 1.3,rmod26= 1.3                      
                                 , sigma = 57499.67)
           , fixed = list( mu = 0) )
