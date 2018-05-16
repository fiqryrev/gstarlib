gstar.w = function(data,round_times=1,jenis=c("IJ","KS","S","B","IJ2","IJ3")){
  #=====================================
  # Checking Conditions
  #=====================================
  if((casefold(jenis)!="ij")&(casefold(jenis)!="ij2")&(casefold(jenis)!="ij3")&(casefold(jenis)!="ks")&(casefold(jenis)!="s")&(casefold(jenis)!="b")) stop(gettextf("weight type '%s' is not supported. Try to use 'IJ' 'IJ2' 'IJ3' 'KS' 'S' or 'B'",jenis),domain=NA)
  if(is.null(data)) stop("Input your data, whether using distance vector or matrix data")
  if(is.data.frame(data)) data = as.matrix.data.frame(data)
  location = dim(data)[2]
  #=====================================
  #Inverse of Distance Weight
  #=====================================
  if(casefold(jenis)=="ij"){
    locs=length(data)
    for(location in 1:locs){
      if(((locs*2)-((location^2)-location))==0) break
      if(location==locs)stop("Dimension of data distance is invalid")
    }
    rMtx = matrix(0,location,location)
    putar = c(0)
    for (i in 1:location){
      for(j in 2:location){
        if((i!=j)&(i<location)&(i<j)){
          putar = putar+1
          rMtx[i,j] = data[putar]
        }
      }
    }
    if(i==location){
      rMtx           = rMtx+t(rMtx)
      rMtxb          = matrix(0,location,location)
      rMtxc          = matrix(0,location,location)
      rMtN           = c()
      for (m in 1:location){
        rMtxb[m,]  = 1/rMtx[m,]
        rMtxb[m,m] = 0
        rMtN[m]    = sum(rMtxb[m,])
        rMtxc[m,]  = rMtxb[m,]/rMtN[m]
        rMtxc[m,m] = 0
      }  
    }
  }  
  #=====================================
  #Square Inverse of Distance Weight
  #=====================================  
  if(casefold(jenis)=="ij2"){
    locs=length(data)
    for(location in 1:locs){
      if(((locs*2)-((location^2)-location))==0) break
      if(location==locs)stop("Dimension of data distance is invalid")
    }
    rMtx = matrix(0,location,location)
    putar = c(0)
    for (i in 1:location){
      for(j in 2:location){
        if((i!=j)&(i<location)&(i<j)){
          putar = putar+1
          rMtx[i,j] = data[putar]
        }
      }
    }
    if(i==location){
      rMtx           = rMtx+t(rMtx)
      rMtxb          = matrix(0,location,location)
      rMtxc          = matrix(0,location,location)
      rMtN           = c()
      for (m in 1:location){
        rMtxb[m,]  = 1/(rMtx[m,]^2)
        rMtxb[m,m] = 0
        rMtN[m]    = sum(rMtxb[m,])
        rMtxc[m,]  = rMtxb[m,]/rMtN[m]
        rMtxc[m,m] = 0
      }  
    }
  }  
  #=====================================
  #Cube Inverse of Distance Weight
  #=====================================    
  if(casefold(jenis)=="ij3"){
    locs=length(data)
    for(location in 1:locs){
      if(((locs*2)-((location^2)-location))==0) break
      if(location==locs)stop("Dimension of data distance is invalid")
    }
    rMtx = matrix(0,location,location)
    putar = c(0)
    for (i in 1:location){
      for(j in 2:location){
        if((i!=j)&(i<location)&(i<j)){
          putar = putar+1
          rMtx[i,j] = data[putar]
        }
      }
    }
    if(i==location){
      rMtx           = rMtx+t(rMtx)
      rMtxb          = matrix(0,location,location)
      rMtxc          = matrix(0,location,location)
      rMtN           = c()
      for (m in 1:location){
        rMtxb[m,]  = 1/(rMtx[m,]^3)
        rMtxb[m,m] = 0
        rMtN[m]    = sum(rMtxb[m,])
        rMtxc[m,]  = rMtxb[m,]/rMtN[m]
        rMtxc[m,m] = 0
      }  
    }
  }  
  #=====================================
  #Binary Weight
  #=====================================    
  if(casefold(jenis)=="b"){
    if(i==location){
      rMtx           = matrix(0,location,location)
      rMtxb          = matrix(0,location,location)
      rMtxc          = matrix(0,location,location)
      rMtN           = c()
      for (m in 1:location){
        rMtxb[m,]  = 1
        rMtxb[m,m] = 0
        rMtN[m]    = sum(rMtxb[m,])
        rMtxc[m,]  = 1
        rMtxc[m,m] = 0
      }  
    }  
  }
  #=====================================
  #Crosscorrelation Weight
  #=====================================    
  if(casefold(jenis)=="ks"){
    if(!is.matrix(data)&!is.data.frame(data)) stop("Cross Correlation should be filled by real data")
    rMtx = matrix(0,location,location)
    for (i in 1:location){
      for(j in 2:location){
        if((i!=j)&(i<location)&(i<j)){
          rMtx[i,j] = ccf(data[,i], data[,j], lag = 1, plot = FALSE)$acf[3]
        }
      }
      for(k in 1:location-1){
        if((i!=k)&(k>=1)&(i>k)){
          rMtx[i,k] = ccf(data[,i], data[,k], lag = 1, plot = FALSE)$acf[3]
        }
      }
      if(i==location){
        rMtxb          = matrix(0,location,location)
        rMtxc          = matrix(0,location,location)
        rMtN           = c()
        for (m in 1:location){
          rMtxb[m,]  = rMtx[m,]
          rMtxb[i,i] = 0
          rMtN[m]    = sum(rMtxb[m,])
          rMtxc[m,]  = rMtxb[m,]/rMtN[m]
          rMtxc[m,m] = 0
        }  
      }
    } 
  }
  #=====================================
  #Uniform Weight
  #=====================================    
  if(casefold(jenis)=="s"){
    rMtxc          = matrix(0,location,location)
    for (m in 1:location){
      rMtxc[m,]  = 1/(location-1)
      rMtxc[m,m] = 0
    }    
  }
  
  #Saving Results
  hasil = round(rMtxc,round_times)
  class(hasil) <- "gstar.w"
  return(hasil)
}