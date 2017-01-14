###############################################################################################
##########################    DATA      CREATION    ###########################################
###############################################################################################


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


#A subset of the Dataset can be seen from the Data below, with time,
#Culmulative Incidence, Culmulative Birth and Culmulative Error.

Fdat = Biweekly(NewData)# Biweekly data
sigmaU = 1
set.seed(143)
U = rnorm(226,0,sd=sigmaU)
Fdat.d = cbind(Fdat[-(227:234),1],Fdat[-(1:8),2],U) # Adjusting for the  delay caused by maternal immunity
# U added was for error, can be ignored
Fdat1 = Cumulative(Fdat.d) # Cumulated table

NewData1 = cbind(Time=time,CIncidence=Fdat1[,1],CBirths=Fdat1[,2], U=Fdat1[,3])

NewData2 = as.data.frame(NewData1 )
head(NewData2)
#error incidence used above if just reports

###############################################################################################
##########################    DATA      CREATED    ############################################
###############################################################################################




###############################################################################################
##########################    MODEL     FITTING    ############################################
###############################################################################################

fit <- lm(CBirths ~ CIncidence, data=NewData2)
summary(fit) # show results



require(locpol)# Package for estimating the parameters of local regression.
d <- data.frame(X = NewData2$CIncidence)
d$Y <- NewData2$CBirths
plot(d$X,d$Y)
sd = sqrt(var(d$X))
h = .3 * sd
lpest1 <- locPolSmootherC(d$X,d$Y , d$X, bw = h , deg = 1, gaussK)
lpest2 <- locpol(Y ~ X, data = d, bw = .3 * sd , kernel = gaussK, deg = 1, xeval = d$X)
lpest2$lpFit
lpest2$residuals

#My codes for Local regression analysis
#######################################
#source('~/GitHub/Measles/locpol.R')
#mylr = lr(d,.3)
#plot(mylr$resd)
#######################################
fitlp = lpest1$beta0 + lpest2$lpFit$Y1 * lpest2$lpFit[,1]
fitlm = fitted.values(fit)
mydat = data.frame(lpest2$lpFit[,1:2], b0=lpest1$beta0, b1=lpest2$lpFit$Y1, 
                   fit=fitlp, resd=lpest2$residuals, fitlm)
###########################################################################

#Data needed is in Fdat.d and mydat

Data = data.frame( Fdat.d[,-3], mydat )
Data[,9] = Data[,1] * Data[,6] #Creating I(t)
Data[,10] = c( NA, Data[,8])[1:226] #Creating Z(t-1)
Data[,11] = c( NA, Data[,9])[1:226] #Creating I(t-1)
colnames(Data)[c(1:2,8,9,10,11)] = c("C","B","Z","I","Zt-1","It-1")


source('~/GitHub/Measles/TSIR parameter estimation.R', echo=F)

Fitit(10) 













