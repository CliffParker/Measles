#Packages
require(pracma)
#Cost function
parcost = function(alpha1 = 1, zeta =1.573e-6, Z0 = -11259.898, I0 =11389.669, rmod = rep(1.349859,26),dta = Data){

  r = rep(NA,226)
  rmod = 1:26
  r[1:226] = rep( rmod , 9)[1:226]
  
  dta[1,10] = Z0
  dta[1,11] = I0
  
  Y = log(dta$I)
  Ycap = log(r) + alpha1 * log(dta$`It-1`) + zeta * dta$`Zt-1`  # without error
  
  ESS = 0
  error = Y - Ycap #assuming out[,2] is for incidence
  ESS = sum(error^2)
  return(ESS)
  
}

#Parameter estimation section
Model_fit=function(par = c( alpha1 = 1, 
                            zeta = 1.573e-6, Z0 = -11259.898, I0 = 11389.669, rmod = rep(1.349859,26))){
  estpars<- fminsearch(parcost, pars <- par, 
            minimize = T, dfree = F,
                       maxiter = 1000, tol = .Machine$double.xmin)
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

