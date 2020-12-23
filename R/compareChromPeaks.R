library(xcms)
source("/home/chuwang/svnrepos/R/msisotope.R")

table1 <- read.table("file.list",header=T)
table2 <- read.table("mass.list",header=T)
table3 <- read.table("probe.list",header=T)

files <- as.list(as.character(table1[,"file"]))
names(files) <- names <- files
charges <- seq(1,6)
masses <- table2[,"mass"]
sequences <- table2[,"sequence"]
probes <- table3[,"probe"]

ppm <- 15*1e-6
Hplus <- 1.0072765

#postscript("output.ps",horizontal=T)
postscript("compareChromPeaks.output.ps",paper="letter")
nrow=ceiling(length(charges)/2)
par(mfrow=c(nrow,length(files)*length(probes)), oma=c(0,0,5,0))
for (i in 1:length(files)) {
  files[[i]] <- xcmsRaw( names(files)[i], profstep=0)
}
for ( m in 1:length(masses) ) { ## per peptide mass
  start.charge <- 1
  mass <- masses[m]
  seq <- sequences[m]
  for ( j in 1:length(charges) ){ ## per charge stages
    charge <- charges[j]
    for ( p in 1:length(probes) ) { ## per probe modifications
      target.mass <- mass + probes[p]
      for ( i in 1:length(files) ) { ## per experiments

        xfile <- files[[i]]
        name <- names[[i]]
        bname <- tail( unlist( strsplit(name,"/") ), 1 )
        mz <- target.mass/charge + Hplus
        mr <- c( mz*(1-ppm), mz*(1+ppm) )
        raw.ECI <- rawEIC(xfile, mr)
        raw.ECI.rt <- xfile@scantime[raw.ECI[[1]]]/60
        ylimit <- max(raw.ECI[[2]])
        plot(raw.ECI.rt, raw.ECI[[2]], type="l",xlab="RT(min)", ylab="intensity")

        mtext(paste("mz: ", as.character(round(mz, digits=5)), "; NL: ",
                    formatC(ylimit, digits=2, format="e")),line=2,cex=0.7)
        mtext(bname,line=1.0,cex=0.7)
        noise <- estimateChromPeakNoise(raw.ECI[[2]])
        if (noise < 10) next
        raw.spec <- cbind(raw.ECI.rt, raw.ECI[[2]])
        dimnames(raw.spec)[[2]] <- c("rt","intensity")
        peaks <- findChromPeaks( raw.spec, noise)
        if (nrow(peaks)>0) {
          for ( k in 1:nrow(peaks) ) {
            rt <- peaks[k,"rt"]
            rt.min <- peaks[k,"rt.min"]
            rt.max <- peaks[k,"rt.max"]
            peak.scan.num <- raw.ECI[[1]][raw.ECI.rt == rt ]
            peak.scan.int <- raw.ECI[[2]][raw.ECI.rt == rt ]
            if ( peak.scan.int < 1e4 ) next
            peak.scan <- getScan(xfile, peak.scan.num,
                                 massrange=c(target.mass/charge-2, target.mass/charge+4 ))
            predicted.dist <- isotope.dist( averagine.count(target.mass) )
            mono.check <- checkChargeAndMonoMass( peak.scan, target.mass, charge, ppm, predicted.dist)
            if (mono.check > 0.8) {
              lines(c(rt,rt),c(ylimit*0.9,ylimit),col="green",lty=1)
            }
          }
        } ## k
      } ## i
    } ## p
    if ( charge%%nrow == 0 ) {
      mtext(paste(m, "Peptide: ", seq, "; Base Mass:", as.character(round(mass,5)),
                  "; Charge:",start.charge, "-", charge), line=2, outer=T)
      start.charge <- charge + 1
    }
  } ## j
} ## m

dev.off()

