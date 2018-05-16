gstar.fit.ols = function(x,arord,dif,W,insample,intercept=FALSE,lbd){
  #================================================================
  #     Fit GSTAR Model with Ordinary Least Square Function       #
  #================================================================
  d_dtx    = dim(x)[1]
  d_lokasi = dim(x)[2]
  #================================================================
  #                      Generating Values                        #
  #================================================================   
  order.min = 1L
  order.max = min(d_dtx - 1L, floor(10 * log10(d_dtx)))
  if(arord>0&arord<=order.max){loopord = arord}else{loopord=order.max}
  X = VarA = VarE = seA = Yhat = A = E = list()
  xaic = xaicc = xbic = xsic = rep.int(Inf, loopord)
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
    if(z==dd1) ZT = c(as.matrix(unlist(ZT)))   #Remake Y
  }
  ztp    = vt = c()
  ZVList = list()
  orde.ar= -1L
  Del1   = TRUE
  for(i in order.min:(order.max+1)){
    orde.ar      = orde.ar+1
    ZVdimension1 = length(ZT)[1]
    ZVdimension2 = (d_lokasi*orde.ar*lbd)
    ztp = data.frame(as.matrix(x[((1+dif):(d1-orde.ar)),]))
    vtp = data.frame(as.matrix(ztp)%*%W)
    d2  = dim(ztp)
    mat1p = matrix(0, prod(dim(ztp)), ncol(ztp))
    mat2p = matrix(0, prod(dim(vtp)), ncol(vtp))
    mat1p[cbind(seq_len(prod(dim(ztp))), rep(seq_along(ztp), each=nrow(ztp)))] = unlist(ztp)
    mat2p[cbind(seq_len(prod(dim(vtp))), rep(seq_along(vtp), each=nrow(vtp)))] = unlist(vtp)
    mat3p = as.matrix(data.frame(mat1p, mat2p))
    
    if(dim(mat3p)[1]==ZVdimension1){
      ZVList[[i]] = mat3p}else{
        if(orde.ar==0){ZVList[[i]]= as.matrix(data.frame("Intercept"=matrix(1,length(ZT))))}else{
          mat3b = matrix(0,(d_lokasi*(orde.ar-1)),dim(mat3p)[2])
          ZVList[[i]]=rbind(mat3b,mat3p)}}
    
    #Checking "null" in ZVList
    checkmate = nullify = c(0L)
    for(ii in order.min:i){
      nullify[ii]  = is.null(ZVList[[ii]])
      checkmate = checkmate+sum(nullify)
      if(ii==i&checkmate>0){
        removenull = ZVList[[which(nullify==1)]] = NULL}}
    
    #Removing First Intercept on non-intercept model
    if (intercept==FALSE&orde.ar==1&Del1==TRUE){ZVList[[1]] = NULL;Del1 = FALSE}
    
    #Set Data Matrix
    if (intercept==TRUE){
      if(orde.ar==0){X[[i]] = as.matrix(data.frame("Intercept"=matrix(1,length(ZT))))}else{
        X[[i]]         = as.matrix(as.data.frame(ZVList))}}
    if(intercept==FALSE){
      if(orde.ar==0){X[[i]] = as.matrix(data.frame("Intercept"=matrix(0,length(ZT))))}else{
        X[[i]]        = as.matrix(as.data.frame(ZVList))}}
    
    XX       = t(X[[i]])%*%X[[i]]
    rank     =  qr(XX)$rank
    if (rank != nrow(XX)) {
      warning(paste("model order: ", i, "singularities in the computation of the projection matrix", 
                    "results are only valid up to model order", i - 1L), domain = NA)
      XXI = XX
    }else{XXI = solve(XX)}
    
    A[[i]]   = ZT%*%X[[i]]%*%XXI
    Yhat[[i]]= A[[i]]%*%t(X[[i]])
    E[[i]]   = (ZT-Yhat[[i]])
    wll=rep.int(1,d1)
    k=2
    VarE[[i]]= as.matrix(tcrossprod(E[[i]])/length(ZT))
    VarA[[i]]= XXI%x%VarE[[i]]
    seA[[i]] = sqrt(diag(VarA[[i]]))
    
    if(intercept==FALSE&orde.ar==0){
      E[[i]] = rep(0,length(ZT))
      loglikhit= NaN
      xaic[i] = NaN 
    }else{
      loglikhit=.5* (sum(log(wll)) - d1 * (log(2 * (22/7)) + 1 - log(d1) +
                                             log(sum(E[[i]]^2))))
      xaic[i] = -2*as.numeric(loglikhit)+k*(ncol(A[[i]])+1)
    }
  }
  IC       = data.frame(AIC=xaic)
  #Remove Inf
  infv     = is.infinite(xaic)
  xaic[which(infv==1)] = NA
  terkecil = arord+1
  #OLS Result  
  X        = X[[terkecil]]
  A        = A[[terkecil]]
  E        = E[[terkecil]]
  Yhat     = matrix(Yhat[[terkecil]],(dim(X)[1]/d_lokasi),d_lokasi)
  VarE     = VarE[[terkecil]]
  VarA     = VarA[[terkecil]]
  seA      = seA[[terkecil]]
  NAIC     = xaic[[terkecil]]
  sdE      = sd(E)
  Emat     = matrix(E,(dim(X)[1]/d_lokasi),d_lokasi)  
  #Create Phi Matrix (GSTAR Coefficients)
  PhijL1 = PhijL2 = M = vector('list',terkecil)
  for(i in 1:arord){
    if(intercept){
      PhijL1[[i]] = diag(A[((i-1)*lbd*d_lokasi+2):((i-1)*lbd*d_lokasi+d_lokasi+1)])
      PhijL2[[i]] = diag(A[((i-1)*lbd*d_lokasi+(d_lokasi+2)):((i-1)*lbd*d_lokasi+lbd*d_lokasi+1)])
    }else{
      PhijL1[[i]] = diag(A[((i-1)*lbd*d_lokasi+1):((i-1)*lbd*d_lokasi+d_lokasi)])
      PhijL2[[i]] = diag(A[((i-1)*lbd*d_lokasi+(d_lokasi+1)):((i-1)*lbd*d_lokasi+lbd*d_lokasi)])
    }  
    M[[i]]      = PhijL1[[i]]+crossprod(PhijL2[[i]],t(W))
  }
  #Saving the Results
  hasil = list()
  hasil$A    = A
  hasil$sdE  = sdE
  hasil$X    = X
  hasil$XXI  = XXI
  hasil$Yhat = Yhat
  hasil$E    = E
  hasil$VarE = VarE
  hasil$VarA = VarA
  hasil$seA  = seA
  hasil$xaic = xaic
  hasil$xaicc= xaicc
  hasil$xbic = xbic
  hasil$xsic = xsic
  hasil$Emat = Emat
  hasil$M    = M
  hasil$IC   = IC
  hasil$sigma= sdE
  return(hasil)
}  