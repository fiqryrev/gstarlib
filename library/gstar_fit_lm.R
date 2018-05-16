gstar.fit.lm = function(x,arord,dif,W,insample,intercept=FALSE,lbd){
  #================================================================
  #       Fit GSTAR Model with Linear Models (Regression)        #
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
  orde.ar= -1L
  Del1   = TRUE
  
  ZVdimension1 = length(ZT)[1]
  ZVdimension2 = (d_lokasi*arord*lbd)
  ztp = data.frame(as.matrix(x[((1+dif):(d1-arord)),]))
  vtp = data.frame(as.matrix(ztp)%*%W)
  d2  = dim(ztp)
  mat1p = matrix(0, prod(dim(ztp)), ncol(ztp))
  mat2p = matrix(0, prod(dim(vtp)), ncol(vtp))
  mat1p[cbind(seq_len(prod(dim(ztp))), rep(seq_along(ztp), each=nrow(ztp)))] = unlist(ztp)
  mat2p[cbind(seq_len(prod(dim(vtp))), rep(seq_along(vtp), each=nrow(vtp)))] = unlist(vtp)
  mat3p = as.matrix(data.frame(mat1p, mat2p))
  X     = mat3p
  for(i in 1:d_lokasi){
    if(i==1)n1=n2=c()
    n1  = c(n1,gettextf("Phi10_Loc%d",i,domain=NA))
    n2  = c(n2,gettextf("Phi11_Loc%d",i,domain=NA))
    if(i == d_lokasi) namae = c(n1,n2)
  }
  colnames(X) = namae
  if(intercept==TRUE){
    lmfit=lm(ZT~X)
  }else{ lmfit=lm(ZT~X-1)}
  A      = lmfit$coefficients
  E      = lmfit$residuals
  Yhat     = matrix(lmfit$fitted.values,(dim(X)[1]/d_lokasi),d_lokasi)
  Emat     = matrix(E,(dim(X)[1]/d_lokasi),d_lokasi)  
  #Create Phi Matrix (GSTAR Coefficients)
  PhijL1 = PhijL2 = M = vector('list',(arord+1))
  for(i in 1:arord){
    PhijL1[[i]] = diag(A[1:d_lokasi])
    PhijL2[[i]] = diag(A[(d_lokasi+1):(lbd*d_lokasi)])
    M[[i]]      = PhijL1[[i]]+crossprod(PhijL2[[i]],t(W))
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