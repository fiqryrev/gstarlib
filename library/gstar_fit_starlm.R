gstar.fit.starlm = function(x,arord,dif,W,insample,intercept=FALSE,lbd){
  #================================================================
  #         Fit STAR Model with Linear Models (Regression)        #
  #================================================================
  d_dtx    = dim(x)[1]
  d_lokasi = dim(x)[2]
  #================================================================
  #                      Generating Values                        #
  #================================================================   
  order.min = 1L
  order.max = min(d_dtx - 1L, floor(10 * log10(d_dtx)))
  if(arord>0&arord<=order.max){loopord = arord}else{loopord=order.max}
  X    = lmfit =Yhat = A = E = list()
  xm   = matrix(0L,1,d_lokasi)
  #================================================================
  #                      Start Analyzing                          #
  #================================================================  
  d1     = dim(x)[1]
  zt     = x[((2+dif):d1),]
  zt1    = x[(1+dif):(d1-1),]
  dd1    = dim(zt)[1]
  ZT = c()
  for (z in 1:dd1){
    ZT = c(ZT,zt[z,])
    if(z==dd1) ZT = c(as.matrix(unlist(ZT)))   #Remake Y/Observations
  }
  ztp    = vt = c()
  ZVList = list()
  
  ZVdimension1 = length(ZT)[1]
  ZVdimension2 = (d_lokasi*arord*lbd)
  ztp = data.frame(as.matrix(x[((1+dif):(d1-arord)),]))
  vtp = data.frame(as.matrix(ztp)%*%W)
  
  ZT1 = VT1 = c()
  for (i in 1:nrow(ztp)){
    ZT1 = c(ZT1,ztp[i,])
    VT1 = c(VT1,vtp[i,])
  }
  
  ZT_Use = as.matrix(c(unlist(ZT1)))
  VT_Use = as.matrix(c(unlist(VT1)))
  X  = cbind("Phi10"=ZT_Use,"Phi11"=VT_Use)
  
  if(intercept==TRUE){
    lmfit = lm(ZT~X)
  }else{
    lmfit=lm(ZT~X-1)
  }
  
  A      = lmfit$coefficients
  E      = lmfit$residuals
  Yhat     = matrix(lmfit$fitted.values,(dim(X)[1]/d_lokasi),d_lokasi)
  Emat     = matrix(E,(dim(X)[1]/d_lokasi),d_lokasi)  
  #Create Phi Matrix (STAR Coefficients)
  PhijL1 = PhijL2 = M = vector('list',(arord+1))
  for(i in 1:arord){
    PhijL1[[i]] = A[arord]
    PhijL2[[i]] = A[arord+1]
    M[[i]]      = PhijL1[[i]]+(PhijL2[[i]]*t(W))
  }
  #Savings the Result
  hasil      = list()
  hasil$A    = A
  hasil$X    = X
  hasil$Yhat = Yhat
  hasil$E    = E
  hasil$Emat = Emat
  hasil$M    = M
  hasil$lmfit= lmfit
  return(hasil)
} 