summary.gstar = function(fit) {
  #Initialization
  takpar=fit$Method
  #Connecting Variables
  if(casefold(takpar)=="lm"|casefold(takpar)=="starlm"){
    lmfit = fit$lmfit
    A    = fit$A
    Yhat = fit$Yhat
    E    = fit$E
    Emat = fit$Emat
    M    = fit$M
  }else{
    A    = fit$A
    sdE  = fit$sdE
    X    = fit$X
    XXI  = fit$XXI
    Yhat = fit$Yhat
    E    = fit$E
    VarE = fit$VarE
    VarA = fit$VarA
    seA  = fit$seA
    xaic = fit$xaic
    Emat = fit$Emat
    M    = fit$M
    IC   = fit$IC
  }
  realdata = fit$realdata
  namalokasi=fit$Namedata
  Dataname = fit$Dataname
  Wname =fit$Nameweight
  Alname=alpha=fit$AlphaUsed
  insample=fit$Usedinsample
  dif=fit$UsedDiff
  W=fit$W
  lbd=fit$Lambda
  intercept=fit$InterceptModel
  takpar=fit$metode
  order=fit$Order
  Orsname=fit$Ordname
  arord = order[1]
  d_data = dim(realdata)[1]
  d_lokasi = dim(realdata)[2]
  realdatain= realdata[(1:insample),]
  realdataout=realdata[(insample+1):(d_data),]
  lbd= fit$Lambda+1
  d1 = dim(realdatain)[1]-dif
  #========================================================================================#
  #                               Parameter Significance                                   #               
  #========================================================================================#
  if(casefold(takpar)=="lm"|casefold(takpar)=="starlm"){
    kesimpulan = summary(lmfit)
  }else{
    #Initial Value
    TValue = Pvalue = Conclusion = Coefs = j = Estimate = Std.Error = Location = c()
    OAR = OL = 0
    OS  = 1 #Spatial Order set only at 1,order 0 and order 1 (lbd=2)
    #Overall Significance Test (F-test)
    RSS    = sum(E^2)
    DF     = length(E)-(lbd*d_lokasi+1L)
    ResVar = RSS/DF
    MSS    = sum(Yhat^2)
    Fhit   = (MSS/(lbd*d_lokasi))/ResVar
    PVF    = df(Fhit,df1=(lbd*d_lokasi),df2=DF)
    if(PVF<alpha){CoF="Significant"}else{CoF="Not Significant"}
    OST    = data.frame("F-Value" = Fhit,"P-Value" = PVF, "Conclusion" = CoF)
    #Partial Significance Test (T-test)
    if(intercept)bykphi=dim(A)[2]-1 else bykphi=dim(A)[2]
    IAR=NAR=NLB=1L;ILB=0L
    for(i in 1:bykphi){
      Coefs[i]=gettextf("Phi%d%d",IAR,ILB,domain=NA)
      Location[i] = namalokasi[NLB]
      NAR = NAR+1;NLB=NLB+1
      if(NLB>d_lokasi){NLB=1;ILB=ILB+1}
      if(ILB>(lbd-1)){ILB=0}
      if(NAR>(lbd*d_lokasi)){NAR=1;IAR=IAR+1}
    }
    bykpar = dim(A)[2]
    for(i in 1:bykpar){
      Estimate[i] = A[i]
      Std.Error[i]= seA[i]
      TValue[i] = A[i]/Std.Error[i]
      Pvalue[i] = pt(abs(TValue[i]),DF,lower.tail=FALSE)
      if(Pvalue[i]<alpha){Conclusion[i]="Significant"}else(Conclusion[i]="Not Significant")
      if(i==(bykpar)){
        if(intercept){
          Location = c("Intercept",Location)
          Coefs    = c(" ",Coefs)}
        hakhir = data.frame(Location,Coefs,Estimate,Std.Error,"T-Value"=TValue,"P-Value"=Pvalue,Conclusion)
      }}
  }
  if(casefold(takpar)=="lm"|casefold(takpar)=="starlm"){
    hasil = list(Parameter=A, Resid = Emat
                 ,Order=Orsname,Lambda=(lbd-1),M=M
                 ,InterceptModel = intercept,Orderan=order, metode=takpar
                 ,Useddata=realdata,Usedweight=W,Usedinsample=insample,UsedDiff=dif
                 ,Namedata = Dataname,Nameweight=Wname,AlphaUsed=Alname)
    cat("**********************************************************************"
        ,"\n*                        Function Space Time                         *"
        ,"\n*                    Autoregressive Model (STAR)                     *"
        ,"\n*                     Lastest Update: 5/16/2018                     *"
        ,"\n**********************************************************************"
        ,"\n","Estimate Method       : ", takpar
        ,"\n**********************************************************************"
        ,"\n","Data to be used       : ", Dataname
        ,"\n","Locations applied     : ", namalokasi
        ,"\n","Weigth to be used     : ", Wname
        ,"\n","Intercept model       : ", intercept
        ,"\n","STAR(p,d,lbd) order   : ", Orsname
        ,"\n","Significance Level    : ", Alname,"%")
    if(dif>1){
      cat("\n","Differencing times    : ", dif, "times")
    }else{
      cat("\n","Differencing times    : ", dif, "time")  
    }
    cat("\n","In sample border      : ", insample)
    cat("\n**********************************************************************"
        ,"\n*                      Model Significance Test                       *"
        ,"\n**********************************************************************")
    print(kesimpulan)  
    cat("\n**********************************************************************"
        ,"\n*                          Output List                               *"
        ,"\n**********************************************************************"
        ,"\n Estimated Parameter    : Parameter"
        ,"\n Residuals              : Resid"
        ,"\n GSTAR Order            : Order"
        ,"\n**********************************************************************"
        ,"\n*                         End of Function                            *"
        ,"\n**********************************************************************","\n")
    
  }else{
    hasil = list(Parameter=A, Resid = Emat,Sigr=hakhir
                 ,Criterion = IC,Order=Orsname,Lambda=(lbd-1),M=M
                 ,InterceptModel = intercept,Orderan=order, metode=takpar
                 ,Useddata=realdata,Usedweight=W,Usedinsample=insample,UsedDiff=dif
                 ,Namedata = Dataname,Nameweight=Wname,AlphaUsed=Alname)
    cat("**********************************************************************"
        ,"\n*                  Function Generalized Space Time                   *"
        ,"\n*                    Autoregressive Model (GSTAR)                    *"
        ,"\n*                     Lastest Update: 15/16/2018                     *"
        ,"\n**********************************************************************"
        ,"\n","Estimate Method       : ", takpar
        ,"\n**********************************************************************"
        ,"\n","Data to be used       : ", Dataname
        ,"\n","Locations applied     : ", namalokasi
        ,"\n","Weigth to be used     : ", Wname
        ,"\n","Intercept model       : ", intercept
        ,"\n","GSTAR(p,d,lbd) order  : ", Orsname
        ,"\n","Significance Level    : ", Alname,"%")
    if(dif>1){
      cat("\n","Differencing times    : ", dif, "times")
    }else{
      cat("\n","Differencing times    : ", dif, "time")  
    }
    cat("\n","In sample border      : ", insample)
    cat("\n**********************************************************************"
        ,"\n*                      Model Significance Test                       *"
        ,"\n**********************************************************************"
        ,"\n                   Overall Fit Test (F-Statistics)                     "        
        ,"\n**********************************************************************"
        ,"\n")
    print(OST) 
    cat("\n**********************************************************************"
        ,"\n                Partial Coefficient Test (T-Statistics)             "
        ,"\n**********************************************************************"
        ,"\n")     
    print(hakhir)
    cat("\n**********************************************************************"
        ,"\n*                  Best Model Information Criterion                  *"
        ,"\n**********************************************************************"
        ,"\n")
    BMIC = as.matrix(IC[arord+1,])
    colnames(BMIC)="AIC"
    rownames(BMIC)=arord
    print(BMIC)
    cat("\n**********************************************************************"
        ,"\n*                          Output List                               *"
        ,"\n**********************************************************************"
        ,"\n Estimated Parameter    : Parameter"
        ,"\n Residuals              : Resid"
        ,"\n Significance Result    : Sigr"
        ,"\n Table of IC            : Criterion"
        ,"\n GSTAR Order            : Order"
        ,"\n**********************************************************************"
        ,"\n*                         End of Function                            *"
        ,"\n**********************************************************************","\n")
  }
  invisible(hasil)
  class(hasil) = "summary.gstar"
  return(hasil)
}