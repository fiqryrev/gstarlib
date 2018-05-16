gstar	= function(x,order=NULL,W=NULL,insample=NULL,alpha=NULL,
                 intercept=TRUE,method=c("OLS","MLE","GLS","LM","STARLM")){
  #================================================================
  #                   Library Initialization                     #
  #================================================================  
  require(TSA)
  require(timeSeries)
  require(starma)
  require(forecast)
  require(tseries)
  require(MTS)
  require(fGarch)
  require(portes)
  require(MVN)
  require(data.table)
  #================================================================
  #                   Function Initialization                     #
  #================================================================
  mulai     = Sys.time()
  m         = match.call()
  realdata  = x
  dimrd     = dim(realdata)
  d_lokasi  = dim(x)[2]
  d_data    = dim(x)[1]
  namalokasi= names(x)
  lambda    = 2       #Still one spatial order
  #================================================================
  #                         Requirements                          #
  #================================================================
  if(is.null(x)){stop("Insert your multivariate data")}else{
    Dataname  = deparse(substitute(x))}
  
  if(NCOL(x)==1){stop("GSTAR is not used to univariate data")}
  
  if(is.null(insample)){
    warning("Insample border is not already set. It is automatically set as 80% from your data length (Yao and Tan, 2001)")
    insample  = floor(0.80*dim(x))[1]
    Insname   = insample}else{
      Insname   = deparse(substitute(insample))
    }
  if(is.null(method)) {
    warning("Estimate method is not already set. It is automatically set to Ordinary Least Square")
    method = "OLS"
    takpar = method}else{
      takpar = method
    }
  if (is.null(W)){
    warning("Weight is not already set. It is automatically set as uniform weight")
    WS = matrix(0,d_lokasi,d_lokasi)
    for (m in 1:d_lokasi){
      WS[m,]  = 1/(d_lokasi-1)
      WS[m,m] = 0
    }
    Wname="Uniform"
    W=WS
  }else{Wname = deparse(substitute(W))} 
  if(is.null(order)|length(order)!=2){
    order = c(1,0)
    dif   = 0L
    arord = 1L
    difname = dif
    Orsname = "GSTAR(1,0,1)"
    warning("GSTAR order is invalid. It is automatically set to GSTAR(1;0;1)")
  }else{
    arord = order[1]
    dif = order[2]
    Difname = dif
    if(casefold(takpar)=="starlm"){
      Orsname = gettextf("STAR(%d,%d,1)",arord,dif,domain=NA) 
    }else{
      Orsname = gettextf("GSTAR(%d,%d,1)",arord,dif,domain=NA)  
    }
  }
  if(is.null(alpha)){
    warning("Alpha is not already set. It is automatically set to 0.05 (5%)")
    alpha     = 0.05
    Alname    = alpha}else{
      Alname    = alpha
    }
  
  if((casefold(method)!="mle")&(casefold(method)!="ols")&(casefold(method)!="gls")&(casefold(method)!="lm")&(casefold(method)!="starlm")) 
    stop(gettextf("Method '%s' is not supported. Try to use 'OLS' 'GLS' 'MLE' 'LM' or 'STARLM'",method,domain=NA))
  
  if (anyNA(x)) stop("There is any NA in your data, try again")
  
  if ((!is.na(insample))&&(isPositiveDefinite(insample))&&(insample<=dimrd[1])){
    x = x[1:insample,]
  }
  
  #================================================================
  #                 Additional Installed Function                 #
  #================================================================  
  
  #Data Centering Function
  pemusatan = function(d){
    blks = dim(d)[2]
    bdt  = dim(d)[1]
    dpus = matrix(0,bdt,blks)
    md = matrix(0,bdt,blks)
    for(i in 1:blks){
      md[,i] = mean(d[,i])
      if(i==blks){
        dpus = sweep(d,2L,md[,i],check.margin = FALSE)
      }
    }
    dpus
  }
  
  det       = function(x) {
    max(0, prod(diag(qr(x)$qr)) * (-1)^(ncol(x) -1))
  }
  
  #Customized Differencing Function
  diff1	    = function(x,o){
    d 	   = dim(x)
    z 	   = data.frame()
    for (i in 1:d[2]){
      z[(1:abs(d[1]-o)),i] = diff(x[,i], differences=o)
    }
    z
  }
  
  #================================================================
  #                     Data Transformation                       #
  #================================================================   
  #  This version is not supporting Automatic Data Transformation #
  #================================================================   
  Transformasi=FALSE
  #================================================================   
  #           Variance Transformation BoxCox, Wei(2006)           #
  #================================================================
  if(Transformasi==TRUE){
    for(i in 1:d_lokasi){
      kpsst = kpss.test(x[,i])
      if(kpsst$p.value>alpha){
        value = BoxCox.lambda(x[,i])  
        if(value!=1L){
          teks  = gettextf("Data %s has been transformed with BoxCox lambda=%e",namalokasi[i],round(value,4),domain=NA)
          x[,i] = BoxCox(x[,i],BoxCox.lambda(x[,i]))
          warning(teks)
        }
      }
    }
  }
  
  #================================================================ 
  #         Mean Transformation by Differencing, Wei(2006)        #
  #================================================================ 
  if ((!is.na(dif))&&(isPositiveDefinite(dif))&&(dif<=d_data)){
    x = diff1(x,dif)
  }
  #================================================================ 
  #                Centering Data for Intercept=TRUE              #
  #================================================================ 
  if(intercept) x = pemusatan(x)
  #================================================================ 
  #                   Parameter Estimation                        #
  #================================================================ 
  if(casefold(takpar)=="mle"){
    fit = gstar.fit.mle(x=x,arord=arord,dif=dif,W=W,insample=insample,intercept=intercept,lbd=lambda)
  }
  if(casefold(takpar)=="ols"){
    fit = gstar.fit.ols(x=x,arord=arord,dif=dif,W=W,insample=insample,intercept=intercept,lbd=lambda)
  }
  if(casefold(takpar)=="gls"){
    fit = gstar.fit.gls(x=x,arord=arord,dif=dif,W=W,insample=insample,intercept=intercept,lbd=lambda)
  }
  if(casefold(takpar)=="lm"){
    fit = gstar.fit.lm(x=x,arord=arord,dif=dif,W=W,insample=insample,intercept=intercept,lbd=lambda)
  }
  if(casefold(takpar)=="starlm"){
    fit = gstar.fit.starlm(x=x,arord=arord,W=W,dif=dif,insample=insample,intercept=intercept,lbd=lambda)
  }
  #================================================================ 
  #                     Getting Parameter                         #
  #================================================================ 
  fit$realdata=realdata
  fit$Dataname=Dataname
  fit$Namedata=namalokasi
  fit$Nameweight=Wname
  fit$AlphaUsed=Alname
  fit$Usedinsample=insample
  fit$UsedDiff=dif
  fit$W=W
  fit$Lambda=(lambda-1)
  fit$InterceptModel=intercept
  fit$metode=takpar
  fit$Order=order
  fit$Ordname=Orsname
  fit$Resid = fit$Emat
  fit$mulai = mulai
  fit$Method = takpar
  fit$arord = arord
  class(fit) <- "gstar"
  invisible(fit)
  #================================================================ 
}