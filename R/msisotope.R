isotope.dist <- function(elements.count, N15.enrichment=1.0) {
  elements <- c( "C", "H", "O", "N", "S", "P", "N15", "H2", "C13", "O18" , "Se", "Cl","Br")
  elements.count.local <- rep(0, length(elements))
  names(elements.count.local) <- elements
  for ( e in names(elements.count)) {
    elements.count.local[e] <- elements.count[e]
  }
  heavy <- c(1.10, 0.015, 0.20, 0.37, 4.21, 0, 100, 100, 100, 100 , 49.16, 24.24, 49.31)/100
  names(heavy) <- elements
  heavy["N15"] <- N15.enrichment

  light <- 1.00 - heavy
  light["Se"] <- 0.2377 # Special for selenium
  #names(elements.count) <- elements
  single.prob <- as.list( elements )
  names(single.prob) <- elements
  all.prob <- numeric(0)
  for ( e in elements ) {
    count <- elements.count.local[e]
    if (count == 0) next
    v <- seq(0,count)
    l <- light[e]
    h <- heavy[e]
    new.prob <- single.prob[[e]] <- round(choose(count,v)*(l^(count-v))*(h^v),4)
    if (e =="O" | e =="S" | e=="Cl" | e=="Br" | e== "O18") { # O, S, Cl and Br isotopes are 2 Da more
      new.prob <- rep(0,2*count+1)
      for( i in 1:(count+1)) {
        new.prob[(i-1)*2+1] <- single.prob[[e]][i]
      }
      single.prob[[e]] <- new.prob
    }
	if (e =="Se") {# Se has the special isotopic distribution
	  single.prob[[e]] <- c(0.0089,0.0000,0.0937,0.0763,0.2377,0.0000,0.4961,0.0000,0.0873)
	  new.prob <- single.prob[[e]]
	}
    #print(single.prob[[e]])
    if ( length(all.prob) > 0 ) {
      d <- length(all.prob)-length(new.prob)
      if ( d > 0) {
        new.prob <- c(new.prob,rep(0,d))
      } else if ( d< 0 ) {
        all.prob <- c( all.prob,rep(0,-d) )
      }
      all.prob <- round(convolve(all.prob, rev(new.prob),type="o"),4)
    } else {
      all.prob <- single.prob[[e]]
    }
  }
  return(all.prob)
##  return(single.prob)
}

## isotope.dist(c(60,13,13,86,2))
averagine.count <- function(input.mass) {
  averagine.mass <- 111.1254
  elements <- c( "C", "H", "O", "N", "S" )
  averagine.comp <- c( 4.9348, 7.7583, 1.4773, 1.3577,0.0417 )
  names(averagine.comp) <- elements
  return( round(averagine.comp*(input.mass/averagine.mass)) )
}

## ms2 triggered?
find.ms2.triggered <- function(xfile, yfile, predicted.mz, rt.range) {
  ms1.scanNums <- which( xfile@scantime>=rt.range[1]&xfile@scantime<=rt.range[2] )
  ms2.matrix <- matrix(0, nrow=length(ms1.scanNums), ncol=length(predicted.mz) )
  dimnames(ms2.matrix)[[1]] <- as.character(ms1.scanNums)
  dimnames(ms2.matrix)[[2]] <- as.character(predicted.mz)
  if (is.null(yfile)) {
    return(ms2.matrix)
  }
  ms2.acq.num.max <- max(xfile@msnAcquisitionNum)
  for (i in 1:dim(ms2.matrix)[1]) {
    ms1.scanNum <- ms1.scanNums[i]
    ms2.acq.num.begin <- xfile@acquisitionNum[ms1.scanNum]+1
    if ( ms1.scanNum < length(xfile@scantime) ) {
      ms2.acq.num.end <- xfile@acquisitionNum[ms1.scanNum+1]-1
    } else {
      ms2.acq.num.end <- ms2.acq.num.max
    }
    ms2.mz <- yfile[yfile[,"num"]>=ms2.acq.num.begin&yfile[,"num"]<=ms2.acq.num.end, "pmz"]
    for (j in 1:dim(ms2.matrix)[2]) {
      ms2.matrix[i,j] <- sum( abs( ms2.mz-predicted.mz[j] )<0.1)
    }
  }
  return(ms2.matrix)
}
##
estimateChromPeakNoise <- function(intensity) {
  mean(intensity)
}

findChromPeaks <- function(spec, noise, sn=5, rtgap=0.2 ) {

##  noise <- estimateChromPeakNoise(spec[,"intensity"] )

  noise <- max(noise, 1e-5)
  spectab <- matrix(nrow = 0, ncol = 4)
  colnames(spectab) <- c("rt", "rt.min", "rt.max", "sn")

  while (spec[i <- which.max(spec[,"intensity"]), "intensity"] > noise*sn) {

    rt <- spec[i,"rt"]
    intensity <- spec[i,"intensity"]
    rt.range <- xcms:::descendValue(spec[,"intensity"], noise, i)
    rt.range <- c(max(1,rt.range[1]-1), min(nrow(spec),rt.range[2]+1) )

    if (rt.range[1] >= 1 && rt.range[2] <= nrow(spec)) {
      rt.min <- spec[rt.range[1],"rt"]
      rt.max <- spec[rt.range[2],"rt"]
      if (!any(abs(spectab[,"rt"] - rt) <= rtgap))
        spectab <- rbind(spectab, c(rt, rt.min, rt.max, spec[i,"intensity"]/noise))
    }
    spec[seq(rt.range[1], rt.range[2]),"intensity"] <- 0
    }

    spectab
}

findPairChromPeaks <- function(rt, light.int, heavy.int, rt.range, local.rt.range, sn=5) {
  #mean of the difference of the ms1 rt 
  rtdiff <- mean( diff(rt) )

  m.light <- cbind(rt,light.int)
  dimnames(m.light)[[2]] <- c("rt","intensity")
  m.light.global <- m.light[rt>=rt.range[1]&rt<=rt.range[2],]
  #noise of the light
  noise.light.global <- estimateChromPeakNoise(m.light.global[,"intensity"])
  m.light.local  <- m.light[rt>=local.rt.range[1]&rt<=local.rt.range[2],]
  noise.light.local <- estimateChromPeakNoise(m.light.local[,"intensity"])
  noise.light <- min( noise.light.global, noise.light.local )
  #find all the peaks of the light
  peaks.light <- findChromPeaks(m.light.global, noise.light, sn, rtgap=0.2)
  #peak num of light
  n.light <- dim(peaks.light)[[1]]

  m.heavy <- cbind(rt,heavy.int)
  dimnames(m.heavy)[[2]] <- c("rt","intensity")
  m.heavy.global <- m.heavy[rt>=rt.range[1]&rt<=rt.range[2],]
  noise.heavy.global <- estimateChromPeakNoise(m.heavy.global[,"intensity"])
  m.heavy.local <- m.heavy[rt>=local.rt.range[1]&rt<=local.rt.range[2],]
  noise.heavy.local <- estimateChromPeakNoise(m.heavy.local[,"intensity"])
  noise.heavy <- min( noise.heavy.global, noise.heavy.local )
  peaks.heavy <- findChromPeaks(m.heavy.global,noise.heavy,sn,rtgap=0.2)
  n.heavy <- dim(peaks.heavy)[[1]]
	#noise of heavy and noise of light
  pair.range <- c(noise.light, noise.heavy)
  if ( n.heavy == 0 | n.light == 0 ) return(pair.range)

  for (i in 1:n.light) {
    rt.i <- peaks.light[i,"rt"]
    d.rt <- abs(rt.i-peaks.heavy[,"rt"])
	#whichi.min returns the position of the minimum number
    j <- which.min(d.rt)
   # if (d.rt[j]>5*rtdiff) next
    rt.j <- peaks.heavy[j,"rt"]
	#for the pair of light and heavy peak, if both of the top points are in the other one's rt.min and rt.max ragion,return the lowest and highest boundary of the pair peaks
    if ( rt.j >=peaks.light[i,"rt.min"] & rt.j <= peaks.light[i,"rt.max"]
        &  rt.i >=peaks.heavy[j,"rt.min"] & rt.i <= peaks.heavy[j,"rt.max"] ) {
      low <- min(peaks.light[i,"rt.min"],peaks.heavy[j,"rt.min"])
	  
      high <- max(peaks.light[i,"rt.max"],peaks.heavy[j,"rt.max"])
      if ( low < high ) {
        pair.range <- c(pair.range,low,high)
      }
    }
  }
  return(pair.range)
}

findSingleChromPeaks <- function(rt, light.int, rt.range, local.rt.range, sn=5) {
  rtdiff <- mean( diff(rt) )

  m.light <- cbind(rt,light.int)
  dimnames(m.light)[[2]] <- c("rt","intensity")
  m.light.global <- m.light[rt>=rt.range[1]&rt<=rt.range[2],]
  noise.light.global <- estimateChromPeakNoise(m.light.global[,"intensity"])
  m.light.local  <- m.light[rt>=local.rt.range[1]&rt<=local.rt.range[2],]
  noise.light.local <- estimateChromPeakNoise(m.light.local[,"intensity"])
  noise.light <- min( noise.light.global, noise.light.local )
  peaks.light <- findChromPeaks(m.light.global, noise.light, sn, rtgap=0.2)
  n.light <- dim(peaks.light)[[1]]

  pair.range <- c(noise.light)
  if ( n.light == 0 ) return(pair.range)

  for (i in 1:n.light) {
    low <- peaks.light[i,"rt.min"]
    high <- peaks.light[i,"rt.max"]
    if ( low < high ) {
      pair.range <- c(pair.range,low,high)
    }
  }
  return(pair.range)
}

readFileFromMsn <- function( xcms.raw ) {
  filename <- xcms.raw@filepath
  filename.base <- strsplit(filename,".mzXML")[[1]][1]
  dtable <- read.table(paste(filename.base,".ms1.mz_int",sep=""), header=T)
  xcms.raw@env$mz <- dtable[,"mz"]
  xcms.raw@env$intensity <- dtable[,"int"]
  dtable <- read.table(paste(filename.base,".ms1.index_table",sep=""), header=T)
  xcms.raw@scanindex <- dtable[,"scanindex"]
  xcms.raw@scantime  <- dtable[,"scantime"]
  xcms.raw@acquisitionNum <- dtable[,"acquisitionNum"]
  if ( length(xcms.raw@msnScanindex) > 0 ) {
    dtable <- read.table(paste(filename.base,".ms2.mz_int",sep=""), header=T)
    xcms.raw@env$msnMz <- dtable[,"mz"]
    xcms.raw@env$msnIntensity <- dtable[,"int"]
    dtable <- read.table(paste(filename.base,".ms2.index_table",sep=""), header=T)
    xcms.raw@msnScanindex <- dtable[, "msnScanindex"]
    xcms.raw@msnAcquisitionNum <- dtable[, "msnAcquisitionNum"]
    xcms.raw@msnPrecursorScan <- dtable[, "msnPrecursorScan"]
    xcms.raw@msnLevel <- dtable[, "msnLevel"]
    xcms.raw@msnRt <- dtable[, "msnRt"]
    xcms.raw@msnPrecursorMz <- dtable[, "msnPrecursorMz"]
    xcms.raw@msnPrecursorIntensity <- dtable[, "msnPrecursorIntensity"]
    xcms.raw@msnPrecursorCharge <- rep(1, length(xcms.raw@msnPrecursorIntensity))
    xcms.raw@msnCollisionEnergy <- rep(xcms.raw@msnCollisionEnergy[1], length(xcms.raw@msnPrecursorIntensity))
  }
  return( xcms.raw )
}

checkChargeAndMonoMass <- function(peak.scan, mono.mass, charge, mz.ppm.cut, predicted.dist) {
  cc <- 0.0
  isotope.mass.unit <- 1.0033548
  Hplus.mass <- 1.0072765
  isomer.max <- which.max(predicted.dist)
  isomer.v <- seq(1, (isomer.max+2) )
  isomer.fit <- predicted.dist[ isomer.v ] / predicted.dist[isomer.max]
  isomer.v <- c(0, isomer.v)
  isomer.fit <- c(0, isomer.fit)
  isomer.mz <- (mono.mass + (isomer.v-1)*isotope.mass.unit)/charge+Hplus.mass
  raw.dist <- isomer.fit
  for ( i in 1:length(raw.dist) ) {
    this.mz <- isomer.mz[i]
    mz.diff <- abs(peak.scan[,1]-this.mz)/this.mz
    if (min(mz.diff) <= mz.ppm.cut ) {
      raw.dist[i] <- peak.scan[which.min(mz.diff),2]
    } else {
      raw.dist[i] <- cc
    }
  }
  max.raw.dist <- max(raw.dist)
  if ( max.raw.dist > 0.0 ) {
    raw.dist <- raw.dist / max.raw.dist
  } else {
    ## no signal found at all
    return(cc)
  }
  ## a strong peak left to the monoisotopic peak
  if (raw.dist[1] > 0.5) return(cc)
  ## only one point for correlation analysis
  npts.expect <- sum(isomer.fit > 0.10)
  if (sum(raw.dist>0) < npts.expect ) return(cc)
  ## quick check on charge states -- take spectrum above 10% intensity cutoff and measure peak gap
  ## not sufficient for overlapping peaks
  i.mz <- which.max(raw.dist)
  if (i.mz == length(raw.dist) ) return(cc)
  mz.range <- c(isomer.mz[i.mz]*(1-mz.ppm.cut), isomer.mz[i.mz+1]*(1+mz.ppm.cut))
  int.range <- max.raw.dist*c(raw.dist[i.mz+1], raw.dist[i.mz])
  tmp1 <- peak.scan[,"intensity"]>=int.range[1] & peak.scan[,"mz"] >= mz.range[1] & peak.scan[,"mz"] <= mz.range[2]
  peak.mz <- peak.scan[tmp1,"mz"]
  n.extra.peak <- length(peak.mz)-2
  wrong.charge <- FALSE
  if (n.extra.peak > 0) {
    for ( j in 1:n.extra.peak ) {
      new.mz <- seq(mz.range[1],mz.range[2],length=j+2)
      wrong.charge <- TRUE
      for ( k in 1:length(new.mz) ) {
        this.mz <- new.mz[k]
        mz.diff <- abs(peak.mz-this.mz)/this.mz
        if (sum(mz.diff<mz.ppm.cut) == 0) {
          wrong.charge <- FALSE
          break
        }
      }
      if (wrong.charge) return(cc)
    }
  }
  cc <- cor(isomer.fit,raw.dist)
  return(cc)
}

multiTitle <- function(...){
###
### multi-coloured title
###
### examples:
###  multiTitle(color="red","Traffic",
###             color="orange"," light ",
###             color="green","signal")
###
### - note triple backslashes needed for embedding quotes:
###
###  multiTitle(color="orange","Hello ",
###             color="red"," \\\"world\\\"!")
###
### Barry Rowlingson <b.rowlingson@lancaster.ac.uk>
###
  l = list(...)
  ic = names(l)=='color'
  colors = unique(unlist(l[ic]))

  for(i in colors){
    color=par()$col.main
    strings=c()
    for(il in 1:length(l)){
      p = l[[il]]
      if(ic[il]){ # if this is a color:
        if(p==i){  # if it's the current color
          current=TRUE
        }else{
          current=FALSE
        }
      }else{ # it's some text
        if(current){
          # set as text
          strings = c(strings,paste('"',p,'"',sep=""))
        }else{
          # set as phantom
          strings = c(strings,paste("phantom(\"",p,"\")",sep=""))
        }
      }
    } # next item
    ## now plot this color
    prod=paste(strings,collapse="*")
    express = paste("expression(",prod,")",sep="")
    e=eval(parse(text=express))
    title(e,col.main=i)
  } # next color
  return()
}

multiMtext <- function(...){
###
### multi-coloured mtext
###
### examples:
###  multiMtext(color="red","Traffic",
###             color="orange"," light ",
###             color="green","signal")
###
### - note triple backslashes needed for embedding quotes:
###
###  multiTitle(color="orange","Hello ",
###             color="red"," \\\"world\\\"!")
###
### Barry Rowlingson <b.rowlingson@lancaster.ac.uk>
###
  l = list(...)
  ic = names(l)=='color'
  colors = unique(unlist(l[ic]))

  for(i in colors){
    color=par()$col.main
    strings=c()
    for(il in 1:length(l)){
      p = l[[il]]
      if(ic[il]){ # if this is a color:
        if(p==i){  # if it's the current color
          current=TRUE
        }else{
          current=FALSE
        }
      }else{ # it's some text
        if(current){
          # set as text
          strings = c(strings,paste('"',p,'"',sep=""))
        }else{
          # set as phantom
          strings = c(strings,paste("phantom(\"",p,"\")",sep=""))
        }
      }
    } # next item
    ## now plot this color
    prod=paste(strings,collapse="*")
    express = paste("expression(",prod,")",sep="")
    e=eval(parse(text=express))
    mtext(e,col=i,line=0.5,outer=T)
  } # next color
  return()
}
