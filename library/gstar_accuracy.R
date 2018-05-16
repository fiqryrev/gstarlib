accuracy.gstar = function(fit){
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
  d_data   = dim(realdata)[1]
  d_lokasi = dim(realdata)[2]
  namalokasi=fit$Namedata
  Dataname = fit$Dataname
  insample=fit$Usedinsample
  dif=fit$UsedDiff
  order=fit$Order
  Orsname=fit$Ordname
  arord = order[1]
  alpha=fit$AlphaUsed
  #====================================
  #Creating Prediction
  #====================================
  realdatain= realdata[(1+dif):(insample),]
  realdataout=realdata[(insample+1):(d_data),]
  rd_matl  = list()
  rd_mat   = as.matrix(realdata)
  rd_mat_hit = c() 
  for(i in 1:(arord)){
    rd_matl[[i]] = rd_mat[(dif+1):(d_data-(i-1)),]
    matNA = matrix(NA,i,d_lokasi)
    if(i>=2)matNA1 = matrix(NA,i-1,d_lokasi)
    if(i>1){
      rd_mat_hit = rd_mat_hit+t(M[[i]]%*%t(rbind(matNA1,rd_matl[[i]])))-
        rbind(matNA,t(M[[i]]%*%t(rd_matl[[i]][-nrow(rd_matl[[i]]),])))  
    }else{
      if(i==1){
        rd_mat_hit = rd_matl[[i]]+t(M[[i]]%*%t(rd_matl[[i]]))-
          rbind(NA,t(M[[i]]%*%t(rd_matl[[i]][-nrow(rd_matl[[i]]),])))
      }}
    
    if(i==arord)rd_mat_hit=data.frame(na.omit(rd_mat_hit))
  }
  Yhatinsample = na.omit(rd_mat_hit[1:(insample+dif),])
  Yhatoutsample = rd_mat_hit[(insample+dif+1):nrow(rd_mat_hit),]
  #====================================
  #Calculate Residual
  #====================================
  residualfull = realdata[(1+arord+dif):(nrow(realdata)),]-rd_mat_hit
  residualin   = data.frame(na.omit(residualfull[1:(insample+dif-arord),]))
  residualout  = residualfull[(insample+dif+1-arord):nrow(rd_mat_hit),]
  #====================================
  #Overall Test
  #====================================
  #Out-sample Test
  #====================================
  MSE = RMSE = MAE = mdAE = MAPE = RMSPE = RMdSPE = MdAPE = c()
  MSE[2] = mean(unlist(residualout)^2)
  RMSE[2] = sqrt(MSE[2])
  MAE[2] = mean(abs(unlist(residualout)))
  mdAE[2] = median(abs(unlist(residualout)))
  pi = (100*unlist(residualout))/unlist(realdataout[1:nrow(residualout),])
  MAPE[2] = mean(abs(pi))
  MdAPE[2] = median(abs(pi))
  RMSPE[2] = sqrt(mean(pi^2))
  RMdSPE[2] = sqrt(median(pi^2)) 
  #====================================
  #In-sample Test
  #====================================
  MSE[1] = mean(unlist(residualin)^2)
  RMSE[1] = sqrt(MSE[1])
  MAE[1] = mean(abs(unlist(residualin)))
  mdAE[1] = median(abs(unlist(residualin)))
  pi = unlist((100*unlist(residualin))/realdatain)
  MAPE[1] = mean(abs(pi))
  MdAPE[1] = median(abs(pi))
  RMSPE[1] = sqrt(mean(pi^2))
  RMdSPE[1] = sqrt(median(pi^2)) 
  
  #====================================
  #Location Test Outsample
  #====================================
  loc_MSE = loc_RMSE = loc_MAE = loc_mdAE =
    loc_MAPE = loc_MdAPE = loc_RMSPE = loc_RMdSPE = c()
  for(i in 1:d_lokasi){
    pi = (100*residualout[,i])/realdataout[1:nrow(residualout),i]
    loc_MSE[i] = mean(residualout[,i]^2)
    loc_RMSE[i] = sqrt(loc_MSE[i])
    loc_MAE[i] = mean(abs(residualout[,i]))
    loc_mdAE[i] = median(abs(residualout[,i]))
    loc_MAPE[i] = mean(abs(pi))
    loc_MdAPE[i] = median(abs(pi))
    loc_RMSPE[i] = sqrt(mean(pi^2))
    loc_RMdSPE[i] = sqrt(median(pi^2)) 
  }
  #====================================
  #Location Test Insample
  #====================================
  inloc_MSE = inloc_RMSE = inloc_MAE = inloc_mdAE =
    inloc_MAPE = inloc_MdAPE = inloc_RMSPE = inloc_RMdSPE = c()
  for(i in 1:d_lokasi){
    pi = (100*residualout[,i])/realdataout[1:nrow(residualout),i]
    inloc_MSE[i] = mean(residualout[,i]^2)
    inloc_RMSE[i] = sqrt(loc_MSE[i])
    inloc_MAE[i] = mean(abs(residualout[,i]))
    inloc_mdAE[i] = median(abs(residualout[,i]))
    inloc_MAPE[i] = mean(abs(pi))
    inloc_MdAPE[i] = median(abs(pi))
    inloc_RMSPE[i] = sqrt(mean(pi^2))
    inloc_RMdSPE[i] = sqrt(median(pi^2)) 
  }
  
  
  
  OverallAccuracy = data.frame(MSE,RMSE,MAE,mdAE,
                               MAPE,MdAPE,RMSPE,RMdSPE)
  rownames(OverallAccuracy) = c("In-sample","Out-sample")
  LocationAccuracy= rbind(loc_MSE,loc_RMSE,loc_MAE,loc_mdAE,
                          loc_MAPE,loc_MdAPE,loc_RMSPE,loc_RMdSPE)
  colnames(LocationAccuracy)=namalokasi
  
  InLocationAccuracy = rbind(inloc_MSE,inloc_RMSE,inloc_MAE,inloc_mdAE,
                             inloc_MAPE,inloc_MdAPE,inloc_RMSPE,inloc_RMdSPE)
  colnames(InLocationAccuracy)=namalokasi
  #XXXXXXXXXXXXXXXXXXXXXXXXXX#
  #     Tracking Signal      #
  #XXXXXXXXXXXXXXXXXXXXXXXXXX#
  if (d_lokasi%%2==0){
    a=d_lokasi/2
    b=d_lokasi/2
  }else{
    a=ceiling(d_lokasi/2)
    b=floor(d_lokasi/2)
  }
  plot.new() #Refresh
  par(mfrow=c(a,b))
  TS = list()
  for(i in 1:d_lokasi){
    y=c(0)
    TS[[i]] = ifelse(cumsum(abs(residualout[,i]-y))==0, 0, 1:length(residualout[,i])*cumsum(residualout[,i] - y)/cumsum(abs(residualout[,i]-y)))
    t = c(1:nrow(residualout))
    plot.ts(y=TS[[i]],x=t,type="l",main=gettextf("Tracking Signal on %s",namalokasi[i],domain=NA),ylim=c(-5,5))
    abline(h=4,col="red",lty=2)   #Batas atas
    abline(h=0,col="black",lty=2) #Titik nol
    abline(h=-4,col="red",lty=2)  #Batas bawah
  }
  TSP = recordPlot() #Saving The Tracking Signal Plot
  
  #XXXXXXXXXXXXXXXXXXXXXXXXXX#
  #Actual vs Prediction Graph#
  #XXXXXXXXXXXXXXXXXXXXXXXXXX#
  if (d_lokasi%%2==0){
    a=d_lokasi/2
    b=d_lokasi/2
  }else{
    a=ceiling(d_lokasi/2)
    b=floor(d_lokasi/2)
  }
  if(d_lokasi==2){
    a=2
    b=1
  }
  plot.new()
  par(mfrow=c(a,b))
  rd_matNA   = matrix(NA,(nrow(realdata)-nrow(rd_mat_hit)),d_lokasi)
  colnames(rd_matNA) = namalokasi
  rd_matplot = rbind(rd_matNA,rd_mat_hit)
  for(i in 1:d_lokasi){
    upper = data.frame(na.omit(rd_matplot[,i])) + (qnorm(1-(alpha/2)))*sigma
    lower = data.frame(na.omit(rd_matplot[,i])) - (qnorm(1-(alpha/2)))*sigma
    plot.ts(realdata[,i], type="n", ylim=range(lower,upper),
            main="Actual vs Prediction",xlab="Series",ylab="Data")
    labels = c("Real Data","Prediction")
    colors = c("black","red")
    legend("bottomright", inset=.005, title=names(realdata)[i],
           labels, lwd=2, col=colors,pch=1,cex=0.8)
    lines(rd_matplot[,i],col='red',lwd=2)
    lines(realdata[,i],lwd=2)
  }
  #Save Results
  APG                     = recordPlot()
  hasil                   = list()
  hasil$Overall_Accuracy  = OverallAccuracy
  hasil$Location_Accuracy = LocationAccuracy
  hasil$Prediction        = rd_mat_hit
  hasil$Resid             = residualfull
  hasil$TSP               = TSP
  hasil$APG               = APG
  hasil$InLocation_Accuracy = InLocationAccuracy
  cat("\n**********************************************************************"
      ,"\n                       Overall Accuracy Check                         "        
      ,"\n**********************************************************************"
      ,"\n")
  print(OverallAccuracy)
  cat("\n**********************************************************************"
      ,"\n                   Out-sample  Locations Accuracy Check               "        
      ,"\n**********************************************************************"
      ,"\n")      
  print(LocationAccuracy)
  cat("\n**********************************************************************"
      ,"\n                 Source: Hyndman and Koehler (2006)                   "        
      ,"\n**********************************************************************"
      ,"\n MSE    = Mean Squared Error                  = mean(e^2)"
      ,"\n RMSE   = Root Mean Squared Error             = sqrt(mean(e^2))"
      ,"\n MAE    = Mean Absolute Error                 = mean(abs(e))"
      ,"\n mdAE   = Median Absolute Error               = median(abs(e))"
      ,"\n MAPE   = Mean Absolute Percentage Error      = mean(e/yt)*100%"
      ,"\n MdAPE  = Median Absolute Percentage Error    = median(e/yt)*100%"
      ,"\n RMSPE  = Root Mean Square Percentage Error   = sqrt(mean(e/yt))*100%"
      ,"\n RMdSPE = Root Median Square Percentage Error = sqrt(median(e/yt))*100%"
      ,"\n")     
  cat("\n**********************************************************************"
      ,"\n*                          Saved Results                             *"
      ,"\n**********************************************************************"      
      ,"\n Prediction                   : Prediction "
      ,"\n Residuals                    : Resid"
      ,"\n Overall Accuracy             : Overall_Accuracy"
      ,"\n In Sample Accuracy/Location  : InLocation_Accuracy"
      ,"\n Out Sample Accuracy/Location : Location_Accuracy"
      ,"\n**********************************************************************"
      ,"\n")
  cat("\n**********************************************************************"
      ,"\n*                          Shown Graphics                            *"
      ,"\n**********************************************************************"      
      ,"\n Actual vs Prediction         : AGP"
      ,"\n Tracking Signal              : TSP"
      ,"\n**********************************************************************"
      ,"\n*                         End of Function                            *"
      ,"\n**********************************************************************"
      ,"\n")
  class(hasil)            = "accuracy.gstar"
  invisible(hasil)
}