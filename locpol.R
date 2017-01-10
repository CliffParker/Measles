# I created this function for local polynonmial regression this should be a good
# Estimate to the real deal as much as possible
lr = function(Data,h){
  #h is the bandwith
  n=dim(Data)[1]
  sd = sqrt(var(Data$X))
  
  Output = data.frame(X = 1, Y = 1, b0 = 1, b1 = 1, fit = 1, resd = 1)
  
  #loop
  
  for(i in 1:n){
    x = Data$X[i]
    a = x - sd*h  #Lower bound
    b = x + sd*h  #Upper bound
    
    #Subseting the data set
    
    Dta = Data[ which( (Data$X >= a) & (Data$X <= b)),]
    row.names(Dta) = NULL
    j =  which(Dta$X == x) # new loop element index
    #performiing linear model
    lmf = lm(Y ~ X, Dta)
    #Fitting output
    place = c( Dta[j,], lmf$coefficients, lmf$fitted.values[j], lmf$residuals[j])
    Output[i,] = place
  }
   
  return(Output)
}

#Data = data.frame(X = 1:20, Y = sin(1:20))
#lr(Data,.3)