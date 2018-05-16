predict.gstar = function(fit,n.ahead=NULL,movingpar=FALSE,CI=FALSE){
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
    sigma = summary(fit$lmfit)$sigma
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
    xaicc= fit$xaicc
    xbic = fit$xbic
    xsic = fit$xsic
    Emat = fit$Emat
    M    = fit$M
    IC   = fit$IC
    sigma= fit$sdE
  }
  realdata = fit$realdata
  namalokasi=fit$Namedata
  Dataname = fit$Dataname
  Wname =fit$Nameweight
  Alname=alpha=fit$AlphaUsed
  insample=fit$Usedinsample
  dif=fit$UsedDiff
  W=fit$W
  lbd=fit$Lambda+1
  intercept=fit$InterceptModel
  takpar=fit$metode
  order=fit$Order
  Orsname=fit$Ordname
  arord = order[1]
  d_data = dim(realdata)[1]
  d_lokasi = dim(realdata)[2]
  realdatain= realdata[(1:insample),]
  realdataout=realdata[(insample+1):(d_data),]
  d1 = dim(realdatain)[1]-dif
  #========================================
  #Checking Conditions
  #========================================
  if(attr(fit,"class")!="gstar") stop("Insert GSTAR fit to this function. Try again")
  
  if(is.null(n.ahead)) {
    n.ahead= 6
    warning("Forecast point is not set yet. It is automatically set to 6 n ahead")}
  if(n.ahead<0) stop("Do not insert negative value. Try again")
  #========================================
  #Creating Forecast
  #========================================
  preddata = realdata
  if(movingpar==FALSE){
    predtable = nullmat =  matrix(NA,n.ahead,d_lokasi)
    rownames(predtable)=c((d_data-dif+2):(d_data-dif+1+n.ahead))
    colnames(predtable) = namalokasi
    rownames(nullmat)=c((d_data-dif+2):(d_data-dif+1+n.ahead))
    colnames(nullmat) = namalokasi
    predasli = list()
    pred_hit = c()
    for(i in 1:n.ahead){
      for(j in 1:arord){
        predasli[[j]] = preddata[(nrow(preddata)-(j)):(nrow(preddata)-(j-1)),]
        if(j==1){
          pred_hit = predasli[[j]][nrow(predasli[[j]]),] + 
            t(M[[j]]%*%t(predasli[[j]][nrow(predasli[[j]]),])) - 
            t(M[[j]]%*%t(predasli[[j]][-nrow(predasli[[j]]),]))
        }else{
          pred_hit = pred_hit + 
            t(M[[j]]%*%t(predasli[[j]][nrow(predasli[[j]]),])) - 
            t(M[[j]]%*%t(predasli[[j]][-nrow(predasli[[j]]),]))
        }
      }
      pred_hit_last = pred_hit
      rownames(pred_hit_last) = nrow(realdata)+i
      preddata = rbind(preddata,pred_hit_last)
      predtable[i,] = unlist(pred_hit_last)
    }
  }else{
    predtable = nullmat =  matrix(NA,n.ahead,d_lokasi)
    if(dif){
      rownames(predtable)=c((d_data-dif+2):(d_data-dif+1+n.ahead))
      colnames(predtable) = namalokasi
      rownames(nullmat)=c((d_data-dif+2):(d_data-dif+1+n.ahead))
      colnames(nullmat) = namalokasi
    }else{
      rownames(predtable)=c((d_data+1):(d_data+n.ahead))
      colnames(predtable) = namalokasi
      rownames(nullmat)=c((d_data+1):(d_data+n.ahead))
      colnames(nullmat) = namalokasi      
    }
    
    predasli = list()
    pred_hit = c()
    for(i in 1:n.ahead){
      for(j in 1:arord){
        predasli[[j]] = preddata[(nrow(preddata)-(j)):(nrow(preddata)-(j-1)),]
        if(j==1){
          pred_hit = predasli[[j]][nrow(predasli[[j]]),] + 
            t(M[[j]]%*%t(predasli[[j]][nrow(predasli[[j]]),])) - 
            t(M[[j]]%*%t(predasli[[j]][-nrow(predasli[[j]]),]))
        }else{
          pred_hit = pred_hit + 
            t(M[[j]]%*%t(predasli[[j]][nrow(predasli[[j]]),])) - 
            t(M[[j]]%*%t(predasli[[j]][-nrow(predasli[[j]]),]))
        }
      }
      pred_hit_last = pred_hit
      rownames(pred_hit_last) = nrow(realdata)+i
      preddata = rbind(preddata,pred_hit_last)
      predtable[i,] = unlist(pred_hit_last)
      gstarmove = gstar(preddata,order=order,W=W,insample=insample,alpha=alpha,
                        intercept=intercept,method=takpar)
      M = gstarmove$M
    }
  }
  #========================================
  #Create Prediction Plot
  #========================================
  if (d_lokasi%%2==0){
    a=d_lokasi/2
    b=d_lokasi/2
  }else{
    a=ceiling(d_lokasi/2)
    b=floor(d_lokasi/2)
  }
  plotdata = rbind(realdata,nullmat)
  par(mfrow=c(a,b))
  for(i in 1:d_lokasi){
    upper = preddata[,i] + (qnorm(1-(alpha/2)))*sigma
    lower = preddata[,i] - (qnorm(1-(alpha/2)))*sigma
    plot(plotdata[,i], type="n", ylim=range(lower,upper),xlab="Series",ylab="Data")
    if(CI==TRUE){
      namaCI = gettextf("CI at (1-%s)",as.character(alpha),domain=NA)
      labels = c("Real Data","Prediction",namaCI)
      colors = c("black","red",rgb(0,0,0.6,0.2))
      polygon(c(time(plotdata[,i]),rev(time(plotdata[,i]))), c(upper,rev(lower)), 
              col=rgb(0,0,0.6,0.2), border=F)
    }else{
      labels = c("Real Data","Prediction")
      colors = c("black","red")  
    }
    legend("bottomright", inset=.005, title=names(realdata)[i],
           labels, lwd=2, col=colors,pch=1,cex=0.8)
    lines(preddata[,i],col='red',lwd=2)
    lines(plotdata[,i],lwd=2)
    out <- (plotdata[,i] < lower | plotdata[,i] > upper)
    points(time(plotdata[,i])[out], plotdata[,i][out], pch=19)
    
  }
  #========================================
  #Save the Result
  #========================================
  hasil = list(Datapred=preddata,RealY=realdata,PredY=predtable
               ,Lambda=lbd,Weight=W
               ,AlphaUsed=alpha,InterceptModel=intercept)
  #========================================
  cat("\014") 
  cat("*************************************************"
      ,"\n*           Function Predict Generalized        *"
      ,"\n*         Space Time Autoregressive (GSTAR)     *"
      ,"\n*            Lastest Update: 12/29/2017         *"
      ,"\n*************************************************"
      ,"\n","Estimate Method       : ", takpar
      ,"\n*************************************************"
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
  cat("\n","Prediction Start      : ", d_data+1
      ,"\n","Prediction End        : ", d_data+n.ahead
      ,"\n","Moving Parameter      : ", movingpar
      ,"\n*************************************************"
      ,"\n*               GSTAR Prediction                *"
      ,"\n*************************************************"
      ,"\n")
  print(predtable)
  cat("*************************************************"
      ,"\n*                  Output List                  *"
      ,"\n*************************************************"
      ,"\n Prediction Data    : Datapred"
      ,"\n*************************************************"
      ,"\n*                End of Function                *"
      ,"\n*************************************************","\n")
  class(hasil) <- "predict.gstar"
  invisible(hasil)
}