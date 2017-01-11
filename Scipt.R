# My script


t = seq(0,9, by = 9/(468))[1:(468)] 
IData = read.table(file.choose())
IData = IData[2:469,2]
BData = read.table(file.choose())
BData = (BData[,2])/13 #Data is recorded quaterly, thus the resulting data is weekly
BData = rep(BData, each = 13)


plot(t,BData, type="l", col="blue")

plot(t,IData, type="l", col="red")



time = seq(0,9, by = 9/(226))[1:(226)] 

#Bi-weekly dataset
NewData = cbind(IData,BData)


Biweekly=function(Data){
  n=nrow(Data)/2
  m=ncol(Data) 
  mat = matrix( ,n,m)
  for (i in 0:n - 1 ){
    mat[i + 1,]= rep(0, m)
    for (j in 1:2) {
      x = (2*i)+j
      mat[i+1,] = c(mat[i+1,]) + c(Data[x,])
    }
  }
  return(mat)
}

Cumulative = function(Data){
  n=nrow(Data)
  m=ncol(Data)
  Dta=matrix( ,n,m)
  Dta[1,] = Data[1,]
  for (i in 2:n){
    Dta[i,] = Dta[i-1,] + Data[i,]
  }
  
  return(Dta)
}


#A subset of the Dataset can be seen from the Data below, with time, Culmulative Incidence, Culmulative Birth and Culmulative Error.

Fdat = Biweekly(NewData)

sigmaU = 1
set.seed(143)
U = rnorm(226,0,sd=sigmaU)
Fdat.d = cbind(Fdat[-(227:234),1],Fdat[-(1:8),2],U)

Fdat1 = Cumulative(Fdat.d) # Adjusting for the  delay caused by maternal immunity

NewData1 = cbind(Time=time,CIncidence=Fdat1[,1],CBirths=Fdat1[,2], U=Fdat1[,3])

NewData2 = as.data.frame(NewData1 )
head(NewData2)

fit <- lm(CBirths ~ CIncidence, data=NewData2)
summary(fit) # show results
plot(predict(fit)~CIncidence, data = NewData2)


confint(fit, level=0.95) # CIs for model parameters 
anova(fit) # anova table 
# covariance matrix for model parameters 

# diagnostic plots 
layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
plot(fit)




#Below are plots of the estimated residuals both with noice and without noise.

Resid.noisefree=as.vector(residuals(fit))# residuals
plot(Resid.noisefree~time, type="o", col="blue")
Resid.noise=as.vector(residuals(fit)) + NewData2$U # residuals with noise
plot(Resid.noise~time, type="o", col="blue")

fit.loess <- loess(CBirths ~ CIncidence, data=NewData2, degree = 1, span = 1.49992)
summary(fit.loess) # show results


Resid.loess=as.vector(fit.loess$residuals)+ NewData2$U# residuals
plot(Resid.loess~time, type="o", col="blue")

fit.loess$kd

require(locpol)# Package for estimating the parameters of local regression.
d <- data.frame(X = NewData2$CIncidence)
d$Y <- NewData2$CBirths
h <- denCVBwSelC(log(d$X), kernel = gaussK)
xeval <- log(d$X)
lpest1 <- locPolSmootherC(log(d$X),log(d$Y) , xeval, bw = h , 1, gaussK)
mean(lpest1$beta1[-c(1,2)])
lpest1
plot(log(d$X),log(d$Y))

attach(lpest1)
dev = beta1 * exp(beta0) * (x)^(beta1 - 1)
finalData = cbind(x, dev, Resid.loess)[-(1:2),]
tail(finalData)



source('~/GitHub/Measles/locpol.R')
 mylr = lr(d,.3)

plot(mylr$resd)

#A Data needed is in Fdat.d and mylr

Data = data.frame( Fdat.d[,-3], mylr )
Data[,9] = Data[,1] * Data[,6]
Data[,10] = c( NA, Data[,8])[1:226]
Data[,11] = c( NA, Data[,9])[1:226]
colnames(Data)[c(1:2,8,9,10,11)] = c("C","B","Z","I","Zt-1","It-1")


source('~/GitHub/Measles/TSIR parameter estimation.R', echo=F)

Fitit(10) 






















