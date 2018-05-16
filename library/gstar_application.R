#=======================================================================================================#
#                                           ###################                                         #
#                                           #        @        #                                         #
#                                           #      @ @        #                                         #
#                                           #        @        #                                         #
#                                           #        @        #                                         #
#                                           #      @@@@@      #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                   INSTALL AND APPLY LIBRARY				                              #
#=======================================================================================================#
lib_ins <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)){install.packages(new.pkg, dependencies = TRUE)}
  sapply(pkg, require, character.only = TRUE)
}
list_package = c("TSA","timeSeries","MTS","fBasics","tseries","starma",
                 "fGarch","portes","MVN","data.table","forecast","tseriesChaos")
lib_ins(list_package)
#=======================================================================================================#
#                                           ###################                                         #
#                                           #       @@@       #                                         #
#                                           #      @  @       #                                         #
#                                           #         @       #                                         #
#                                           #        @        #                                         #
#                                           #      @@@@@      #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                          LIBRARY R		                                            #
#=======================================================================================================#
#Why should we use those libraries?
#---------------------------------------------------------
#Library name           Functions
#---------------------------------------------------------
#library(TSA)           #Stationary Test in R
#library(timeSeries)    #Data Timeseries
#library(tseriesChaos)  #Chaos Time Series
#library(starma)        #Stacf & stpacf
#library(MTS)           #Multivariate ARCH Test
#library(fGarch)        #ARCH-GARCH model fit
#library(portes)        #Multivariate White Noise Test
#library(MVN)           #Multivariate Normality Test
#library(forecast)      #Standard TimeSeries 
#library(tseries)       #Standard TimeSeries 
#library(data.table)    #Make Data Table
#library(fBasics)       #RainbowPallete
#---------------------------------------------------------
#=======================================================================================================#
#                                           ###################                                         #
#                                           #       @@@@@     #                                         #
#                                           #          @@     #                                         #
#                                           #       @@@@@     #                                         #
#                                           #          @@     #                                         #
#                                           #       @@@@@     #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                 Input Data	& Summary Statistics		                              #
#=======================================================================================================#
setwd("~sompo_labeling/")
dt  = read.csv("DataIHKJabarJanuari8.csv",header=T)

insample=40
summary(dt)
attach(dt)
banyak_data   = dim(dt)[1]
banyak_lokasi = dim(dt)[2]
datain 	     = dt[1:insample,]
dataout 	   = dt[(insample+1):nrow(dt),]
#Grafik Data Keseluruhan
ts.plot(dt,ylab="",main="Data IHK Provinsi Jawa Barat", col=c("red","green","blue","magenta"),lwd=2)
legend("bottomright",c("Bekasi","Bogor","Sukabumi","Depok"),cex=1, lty = 2, text.font=1, col=c("red","green","blue","magenta"))
names(dt)
#Grafik Data Perlokasi
par(mfrow=c(2,2))			
for(i in 1:banyak_lokasi){
  color = rgb(red=sample(1:255,1),green=sample(1:255,1),blue=sample(1:255,1),max=255)
  teks = c("Data IHK Kota: ",names(dt)[i])
  plot.ts(dt[,i],ylab="",main=teks,col=color,lwd=2)
  abline(h=mean(dt[,i]),col="red",lwd=2)
}
summary(dt);cor(dt)     
#Range	Interquantile Range	Varians	Standar Deviasi	Kurtosis	Skewness
boxplot(dt,col=rainbowPalette(n=5),main="Box Plot IHK di Empat Kota Provinsi Jawa Barat")

#=======================================================================================================#
#                                           ###################                                         #
#                                           #        @ @      #                                         #
#                                           #       @  @      #                                         #
#                                           #      @@@@@@     #                                         #
#                                           #          @      #                                         #
#                                           #          @      #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                    Pemusatan Data in sample 		                                  #
#=======================================================================================================#
attach(datain)
datainpus = matrix(0,dim(datain)[1],banyak_lokasi)
for(i in 1:banyak_lokasi){
  datainpus[,i] = datain[,i]-mean(datain[,i])
  color = rgb(red=sample(1:255,1),green=sample(1:255,1),blue=sample(1:255,1),max=255)
  teks = c("Data IHK Kota: ",names(dt)[i])
  plot.ts(datainpus[,i],ylab="",main=teks,col=color,lwd=2)
  abline(h=mean(datainpus[,i]),col="red",lwd=2)
}
attach(as.data.frame(datainpus))

dV = list()
for(i in 1:banyak_lokasi){
  dV[[i]] = datainpus[,i]
}

#=======================================================================================================#
#                                           ###################                                         #
#                                           #      @@@@@@     #                                         #
#                                           #      @          #                                         #
#                                           #      @@@@@      #                                         #
#                                           #           @     #                                         #
#                                           #      @@@@@      #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                Uji Stasioneritas Data in sample	                                #
#=======================================================================================================#
#Staionary in Mean

#Stationary in Variance within Trend Time Series Model
#KPSS test unit root for data time series with trend
kpss=c()
for (i in 1:banyak_lokasi){
  kpss[[i]] = kpss.test(datain[,i])
}
kpss

#=======================================================================================================#
#                                           ###################                                         #
#                                           #       @@@@      #                                         #
#                                           #      @          #                                         #
#                                           #     @@@@@@      #                                         #
#                                           #      @   @      #                                         #
#                                           #       @@@@      #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                           Differencing Data in sample dan Grafiknya	                            # 
#=======================================================================================================#
dDatainpus = c()
par(mfrow=c(2,2))		
for(i in 1:banyak_lokasi){
  dV[[i]] = diff(dV[[i]])
  color = rgb(red=sample(1:255,1),green=sample(1:255,1),blue=sample(1:255,1),max=255)
  teks = c("Data IHK Kota: ",names(dt)[i])
  plot.ts(dV[[i]],ylab="",main=teks,col=color,lwd=2)
  abline(h=mean(dV[[i]]),col="red",lwd=2)
  dDatainpus = cbind(dDatainpus,dV[[i]])
}
for(i in 1:banyak_lokasi){
  color = rgb(red=sample(1:255,1),green=sample(1:255,1),blue=sample(1:255,1),max=255)
  teks = c("Data IHK Kota: ",names(dt)[i])
  plot.ts(dDatainpus[,i],ylab="",main=teks,col=color,lwd=2)
  abline(h=mean(dDatainpus[,i]),col="red",lwd=2)
}
#=======================================================================================================#
#                                           ###################                                         #
#                                           #      @@@@@@     #                                         #
#                                           #          @      #                                         #
#                                           #         @       #                                         #
#                                           #        @        #                                         #
#                                           #       @         #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                           Uji Stasioneritas Data in sample Differencing	                        # 
#=======================================================================================================#
#Stationary in Variance
kpss = c()
for (i in 1:banyak_lokasi){
  kpss[[i]] = kpss.test(dDatainpus[,i],null="Level")
}
kpss
#=======================================================================================================#
#                                           ###################                                         #
#                                           #      @@@@@@     #                                         #
#                                           #      @    @     #                                         #
#                                           #      @@@@@@     #                                         #
#                                           #      @    @     #                                         #
#                                           #      @@@@@@     #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                     Pembobot WIJ dan WKS	                                        # 
#=======================================================================================================#
inputjarak = c(57.9,112.9,39.3,65.0,43.1,100.1)
WIJ = gstar.w(inputjarak,pembulatan=2,jenis="IJ")  
WIJ
WKS = gstar.w(datainpus,pembulatan=4,jenis="ks") 
WKS
#=======================================================================================================#
#                                           ###################                                         #
#                                           #      @@@@@@     #                                         #
#                                           #      @    @     #                                         #
#                                           #      @@@@@@     #                                         #
#                                           #           @     #                                         #
#                                           #      @@@@@@     #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                     Stpacf Data Stasioner	                                      # 
#=======================================================================================================#
#Pembobot WKS#
stpacf(dDatainpus,WKS, tlag.max=10, plot=c(TRUE), use.ggplot=TRUE)#plot
stpacf(dDatainpus,WKS, tlag.max=10, plot=c(FALSE), use.ggplot=TRUE)#nilai Plot
#=======================================================================================================#
#                                           ###################                                         #
#                                           #   @@     @@@@@@ #                                         #
#                                           #  @ @     @    @ #                                         #
#                                           #    @     @    @ #                                         #
#                                           #    @     @    @ #                                         #
#                                           #  @@@@@   @@@@@@ #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                         SYNTAX GSTAR	                                           # 
#=======================================================================================================#
#Skenario 1
intercept = FALSE;insample=40;alpha=0.05;order=c(1,1)
fitku1 = gstar(x=dt,order=c(1,1),W=WIJ,alpha=0.05,intercept=FALSE,
               insample=insample,method="OLS")
fitku4 = gstar(x=dt,order=c(4,1),W=WIJ,alpha=0.05,intercept=FALSE,
               insample=insample,method="OLS")
fitku5 = gstar(x=dt,order=c(5,1),W=WKS,alpha=0.05,intercept=FALSE,
               insample=insample,method="OLS")
fitku7 = gstar(x=dt,order=c(7,1),W=WKS,alpha=0.05,intercept=FALSE,
               insample=insample,method="OLS")

#Hasil
summary(fitku1)
summary(fitku4)
summary(fitku5)
summary(fitku7)
ac1 = accuracy(fitku1)
ac4 = accuracy(fitku4)
ac5 = accuracy(fitku5)
ac7 = accuracy(fitku7)
#=======================================================================================================#
#                                           ###################                                         #
#                                           #    @@     @@    #                                         #
#                                           #   @ @    @ @    #                                         #
#                                           #     @      @    #                                         #
#                                           #     @      @    #                                         #
#                                           #   @@@@@  @@@@@  #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                               Menghitung Kekeliruan model	GSTAR 	                              # 
#=======================================================================================================#
gks_vtotal	   = as.matrix(fitku1$Resid)

ed_ks          = dim(gks_vtotal)
edt_ks         = ed_ks[1]*ed_ks[2]
gks_v          = matrix(0,ed_ks[1],ed_ks[2])
gks_vk         = matrix(0,ed_ks[1],ed_ks[2])
gks_vkx        = matrix(0,ed_ks[1],ed_ks[2])
for (i in 1:ed_ks[2]){
  gks_v[,i] = as.vector(gks_vtotal[,i])
  gks_vk[,i] = gks_v[,i]^2
  gks_vkx[,i] = gks_vk[,i]
}
cor(gks_vtotal) #melihat efek korelasi
#=======================================================================================================#
#                                           ###################                                         #
#                                           #    @@     @@@   #                                         #
#                                           #   @ @    @   @  #                                         #
#                                           #     @        @  #                                         #
#                                           #     @       @   #                                         #
#                                           #   @@@@@  @@@@@@ #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                            Uji diagnostik model GSTAR                                      # 
#=======================================================================================================#
#Pengujian Homoskedastisitas 
MarchTest(gks_vkx)           
#Pengujian Multivariat Normal
roystonTest(gks_v,qqplot=T) 
#Pengujian White Noise
portest(gks_vkx,lags=seq_len(30),test=c("BoxPierce"), MonteCarlo=FALSE)
#=======================================================================================================#
#                                           ###################                                         #
#                                           #    @@    @@@@   #                                         #
#                                           #   @ @       @   #                                         #
#                                           #     @    @@@@   #                                         #
#                                           #     @       @   #                                         #
#                                           #   @@@@@  @@@@   #                                         #
#                                           ###################                                         #
#=======================================================================================================#
#  		                                        Peramalan GSTAR  	                                        # 
#=======================================================================================================#
predict(fitku1,n.ahead=6,movingpar=TRUE,CI=TRUE)
predict(fitku4,n.ahead=6,movingpar=TRUE,CI=TRUE)
predict(fitku5,n.ahead=6,movingpar=TRUE,CI=TRUE)
predict(fitku7,n.ahead=6,movingpar=TRUE,CI=TRUE)

