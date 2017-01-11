
parset = function(alpha1 = 1, zeta =1.573e-6, Z0 = 5005.874, I0 =11673.423, rmod = rep(.2,26),dta = Data){

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


require(pracma)
Model_fit=function(){
  estpars<- fminsearch(parset, pars <- c( alpha1 = 1, 
            zeta = 1.573e-6, Z0 =5005.874, I0 = 11673.423, rmod = rep(.2,26)), 
            minimize = TRUE, dfree = TRUE,
                       maxiter = 1000, tol = .Machine$double.eps^(2/3))
  return(estpars)
}

Model_fit()
