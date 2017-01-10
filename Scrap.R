install.packages("locpol")
require(locpol)

Dta = data.frame(x = 1:100 , y = sin(pi*(1:100)))

attach(Dta)
Obj = locpol(y~x,Dta,weig=rep(1,nrow(Dta)),bw= 1,kernel=gaussK,deg=2,
             xeval=x)
Obj$lpFit
Obj$residuals