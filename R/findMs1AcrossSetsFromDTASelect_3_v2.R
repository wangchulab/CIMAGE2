library(xcms)
cimage.path <- Sys.getenv("CIMAGE_PATH")
source(paste(cimage.path,"/R/inputparams.R",sep=""))
source(paste(cimage.path,"/R/msisotope.R",sep=""))

## file name from input args
args <- commandArgs(trailingOnly=T)
param.file <- args[1]
args <- args[-1]

## read parameters from a file
params <- read.input.params(param.file)
quantification_type <- c("medium","heavy")

## initialize atom mass table
atom.mass.vec <- init.atom.mass()
## initialize chemical composition table
light.chem.table <- read.chem.table(params[["light.chem.table"]])
medium.chem.table <- read.chem.table(params[["medium.chem.table"]])
heavy.chem.table <- read.chem.table(params[["heavy.chem.table"]])
## initialize amino acid mass table
light.aa.mass <- init.aa.mass(atom.mass.vec, light.chem.table)
medium.aa.mass <- init.aa.mass(atom.mass.vec, medium.chem.table)
heavy.aa.mass <- init.aa.mass(atom.mass.vec, heavy.chem.table)

## ppm to extract chromatographic peaks
mz.ppm.cut <- as.numeric(params[["ppm.tolerance"]]) * 1E-6
## N15 enrichment ratio ##
N15.enrichment <- as.numeric(params[["N15.enrichment"]])
## nature mass difference between C12 and C13
isotope.mass.unit <- atom.mass.vec["C13"] - atom.mass.vec["C"]
## natural mass difference between N14 and N15
isotope.mass.unit.N15 <- atom.mass.vec["N15"] - atom.mass.vec["N"]
# mass of a proton
Hplus.mass <- atom.mass.vec["Hplus"]

## output folder
output.path <- params[["output.path"]]
dir.create(output.path)
## the table with protein names
ipi.name.table <- read.table("ipi_name.table",sep="\t",header=T,comment.char="")
## the table with mass and scan number from DTASelect
cross.table <- read.table("cross_scan.table", header=T, check.names=F,comment.char="")
#cross.table[,"mass"] <- cross.table[,"mass"] + probe.mass
split.table <- matrix(unlist(strsplit(as.character(cross.table[,"key"]),":")), byrow=T,ncol=4)
dimnames(split.table)[[2]] <- c("ipi","peptide","charge","segment")
cross.table <- cbind(cross.table, split.table)
uniq.ipi.peptides <- as.factor(paste(cross.table[,"ipi"], cross.table[,"peptide"],sep=":"))
entry.levels <- levels( uniq.ipi.peptides )
## all_scan.table
all.scan.table <- read.table("all_scan.table", header=T, as.is=T,comment.char="",colClasses=c("character","character","numeric","character"))
## file name tags
cross.vec <- as.character(args)
print(args)
ncross <- length(cross.vec)
print(ncross)
# handle switching of Heavy to light ratio, by default it is light vs heavy
HL.ratios <- rep(FALSE,ncross)
j <- 0
# by default is to cacullate, L/H ratio, if putting the HL behind the file name, the HL.ratios will be true and will calculate the H/L ratio instead
for ( arg in as.character(args) ) {
  j <- j+1
  cross.vec[j] <- sub("_HL$","",arg)
  if ( cross.vec[j] != arg ) {
    HL.ratios[j] <- TRUE
  }
}
## find all matched mzXML input files in upper directory
if(TRUE){
#if(FALSE){
input.path <- getwd()
mzXML.names <- list.files(path="../",pattern="mzML$")
mzXML.files <- as.list( mzXML.names )
names(mzXML.files) <- mzXML.names
for (name in mzXML.names) {
  cat(paste(name,"\n",sep=""))
  mzXML.files[[name]] <- xcmsRaw( paste("../",name,sep=""), profstep=0, includeMSn=T)
}
}
#print(mzXML.files)
## more parameters from input files
## retention time window for alignment across multiple samples
rt.window <- as.numeric(params[["rt.window"]])
rt.window.width <- rt.window * 60
## local retention time window for noise line calculation
local.rt.window <- as.numeric(params[["local.rt.window"]])
local.rt.window.width <- local.rt.window * 60
## signal/noise ratio for peak picking
sn <- as.numeric(params[["sn"]])
### range cutoff for calculated ratios###
ratio.range <- as.numeric(params[["ratio.range"]])
### isotope envelope score cutoff ###
env.score.cutoff <- as.numeric(params[["env.score.cutoff"]])
### coelution profile r2 cutoff ###
r2.cutoff <- as.numeric(params[["r2.cutoff"]])
### minimum peak width in numbers of time points###
minimum.peak.points <- as.numeric(params[["minimum.peak.points"]])
### choose peak pairs with MS2 data only ###
peaks.with.ms2.only <- as.logical(params[["peaks.with.ms2.only"]])
### singleton case ratio ###
singleton.ratio <- as.numeric(params[["singleton.ratio"]])
## column names for calculated ratios
integrated.area.ratio.LH <- paste("IR_LH",cross.vec,sep=".")
integrated.area.ratio.LM <- paste("IR_LM",cross.vec,sep=".")
integrated.area.ratio.HM <- paste("IR_HM",cross.vec,sep=".")
linear.regression.ratio.LH <- paste("LR_LH",cross.vec,sep = ".")
linear.regression.ratio.LM <- paste("LR_LM",cross.vec,sep = ".")
linear.regression.ratio.HM <- paste("LR_HM",cross.vec,sep = ".")
peak.noise.information <- paste("NP",cross.vec,sep=".")
linear.regression.R2.LH <- paste("R2_LH",cross.vec,sep=".")
linear.regression.R2.LM <- paste("R2_LM",cross.vec,sep=".")
linear.regression.R2.HM <- paste("R2_HM",cross.vec,sep=".")
light.integrated.area <- paste("INT",cross.vec,sep=".")
column.names <- c("index","ipi", "description", "symbol", "sequence", "mass", "charge", "segment",
                  integrated.area.ratio.LH,integrated.area.ratio.LM,integrated.area.ratio.HM, linear.regression.ratio.LH, linear.regression.ratio.LM,linear.regression.ratio.HM, 
				  light.integrated.area, peak.noise.information, linear.regression.R2.LH,linear.regression.R2.LM,linear.regression.R2.HM, "entry", "link" )
out.df <- matrix(nrow=0, ncol=length(column.names))
colnames(out.df) <- column.names

## output name
out.filename.base <- paste("output_rt_",as.character(rt.window),"_sn_",
                      as.character(sn),sep="")
#out.filename <- paste(output.path,"/",out.filename.base,".pdf",sep="")
## output layout
#pdf( out.filename, height=4*ncross, width=11, paper="special")
layout.vec <- row.layout.vec <- c(1,1,2,1,1,3,4,4,5,4,4,6,7,7,8,7,7,9)
if ( ncross > 1 ) {
  for (i in 1:(ncross-1)) {
    layout.vec <- c(layout.vec,(row.layout.vec+i*3))
  }
}
#creat the layout matrix for picture drawing
layout.matrix <- matrix(layout.vec,byrow=T,ncol=3)
layout(layout.matrix)
#oma is used to set the boundary of the picture, las is used to set the scale interval
par(oma=c(0,0,5,0), las=0)

dir.create(paste(output.path,"/PNG",sep=""))
npages <- dim(cross.table)[1]
message(paste("Total number of pages are ",npages,sep=""))
#npages <- 11
#for ( i in 133:133) { 
for ( i in 1:npages) {
  #print(i)
  i.folder <- floor((i-1)/500)
  i.page <- (i-1)%%500
  if (! i.page) {
    dir.create(paste(output.path,"/PNG/",i.folder,sep=""))
    message(paste("working on pages ",i.folder*500+1,"--",min((i.folder+1)*500,npages),sep=""))
  }
  out.filename <- paste(output.path,"/PNG/",i.folder,"/",out.filename.base,".",i.folder,"_",i.page,".png",sep="")
  png( out.filename, height=1200*ncross, width=350*3,pointsize=16)
  #out.filename <- paste(output.path,"/",out.filename.base,"_",i,".pdf",sep="")
  #pdf( out.filename, height=4*ncross, width=11, paper="special")
  layout(layout.matrix)
  par(oma=c(0,0,5,0), las=0)
  key <- cross.table[i,"key"]
  tmp.vec <- unlist( strsplit(as.character(key),":") )
  ipi <- tmp.vec[1]
  peptide <- tmp.vec[2]
  entry.tag <- paste( ipi, peptide, sep=":")
  entry.index <- which( entry.tag == entry.levels )
#  if (! grepl('\\*',peptide) & ! grepl('SILAC',param.file)) {
#	peptide.list <- unlist(strsplit(peptide,split="\\."))
#	peptide.list[2] <- paste(peptide.list[2],"*",sep="")
#	peptide <- paste(peptide.list,collapse=".")
#	}
  charge <- as.integer(tmp.vec[3])
  segment <- tmp.vec[4]
  
  
  description <- as.character(ipi.name.table[ipi,"name"])
  symbol<- strsplit(description, " ")[[1]][1]
  ## momo.mass (light and heavy) and mass of most abundant isotopes (light and heavy)
  mono.mass <- calc.peptide.mass( peptide, light.aa.mass)
  #predicted.dist <- isotope.dist( averagine.count(mono.mass) )
  elements.count <- calc.num.elements(peptide, light.chem.table)
  predicted.dist <- isotope.dist(elements.count)
  i.max <- which.max(predicted.dist)
  mass <- mono.mass + (i.max - 1)*isotope.mass.unit
  mono.mass.medium <- calc.peptide.mass( peptide, medium.aa.mass)
  mono.mass.heavy <- calc.peptide.mass( peptide, heavy.aa.mass)
  mass.medium <- mono.mass.medium + mass - mono.mass
  mass.heavy <- mono.mass.heavy + mass - mono.mass
  elements.count.medium <- calc.num.elements(peptide, medium.chem.table)
  predicted.dist.medium <- isotope.dist(elements.count.medium,N15.enrichment)
  elements.count.heavy <- calc.num.elements(peptide, heavy.chem.table)
  predicted.dist.heavy <- isotope.dist(elements.count.heavy,N15.enrichment)
  ## mass delta between light and heavy
  mass_shift.medium <- sum((elements.count.medium-elements.count)[c("N15","H2","C13")])
  correction.factor.medium <- predicted.dist[i.max]/predicted.dist.medium[i.max+mass_shift.medium]
  mass_shift.heavy <- sum((elements.count.heavy-elements.count)[c("N15","H2","C13")])
  correction.factor.heavy <- predicted.dist[i.max]/predicted.dist.heavy[i.max+mass_shift.heavy]
  ## mz
  mono.mz <- mono.mass/charge + Hplus.mass
  mz.light <- mass/charge + Hplus.mass
  mono.mz.medium <- mono.mass.medium/charge + Hplus.mass
  mz.medium <- mass.medium/charge + Hplus.mass
  mono.mz.heavy <- mono.mass.heavy/charge + Hplus.mass
  mz.heavy <-  mass.heavy/charge + Hplus.mass
  ## scan number
  raw.scan.num <- cross.table[i,cross.vec]
  ms1.scan.rt <- ms1.scan.num <- exist.index <- which( raw.scan.num > 0 )
  # do not know what it is for
  for ( k in 1:length(exist.index) ) {
    kk <- exist.index[k]
    raw.file <- paste( cross.vec[kk], "_", segment,".mzML",sep="")
    xfile <- mzXML.files[[raw.file]]
    ms1.scan.num[k] <- which(xfile@acquisitionNum > as.integer(raw.scan.num[kk]))[1]-1
    if (is.na(ms1.scan.num[k])) {
      ms1.scan.num[k] <- length(xfile@acquisitionNum)
    }
    ms1.scan.rt[k] <- xfile@scantime[ms1.scan.num[k]]
  }

  r2.v_HM <- r2.v_LH <- r2.v_LM <- l_ratios_HM <- l_ratios_LH <- l_ratios_LM <- NP.value <- light.int.v <- i_ratios_HM <- i_ratios_LH <- i_ratios_LM <- rep(NA,ncross)
  #print(ncross)
  for ( j in 1:ncross ) {
    #print(mzXML.files)
    raw.file <- paste( cross.vec[j], "_", segment,".mzML",sep="")
	#print(raw.file)
    xfile <- mzXML.files[[raw.file]]
	#print(xfile)
    ## tag * and tag rt line
    if ( j %in% exist.index ) {
      tag <- "*"
      tag_ms1_scan_num <- ms1.scan.num[match(j,exist.index)]
	  #scantime of the ms1 scan in minite
      tag.rt <- xfile@scantime[tag_ms1_scan_num]/60
    } else {
      tag <- ""
      tag_ms1_scan_num <- NA
      tag.rt <- NA
    }
	#print(xfile)
    ##chromatogram bottom; EIC is the  chromatogram of ion (extracted ion chromatogram)
    raw.ECI.light <- rawEIC(xfile, c(mz.light*(1-mz.ppm.cut), mz.light*(1+mz.ppm.cut)) )
	raw.ECI.medium <- rawEIC(xfile, c(mz.medium*(1-mz.ppm.cut), mz.medium*(1+mz.ppm.cut)) )
    raw.ECI.heavy <- rawEIC(xfile, c(mz.heavy*(1-mz.ppm.cut), mz.heavy*(1+mz.ppm.cut)) )
    scan.time.range <- range(xfile@scantime)
	#calculate the left boundary of the rt.window
    rt.min <- min(ms1.scan.rt)-rt.window.width
    if (rt.min > scan.time.range[2]) {
      rt.min <- scan.time.range[2] - 2*rt.window.width
    } else {
      rt.min <- max(rt.min, scan.time.range[1] )
    }
	#calculate the right boundary of the rt.window
    rt.max <- max(ms1.scan.rt)+rt.window.width
    if (rt.max < scan.time.range[1] ) {
      rt.max <- scan.time.range[1] + 2*rt.window.width
    } else {
      rt.max <- min(rt.max,scan.time.range[2])
    }
    if ( (rt.max - rt.min) < 2*rt.window.width ) {
      if ( rt.max == scan.time.range[2] ) {
        rt.min <- rt.max - 2*rt.window.width
      } else if (rt.min == scan.time.range[1]) {
        rt.max <- rt.min + 2*rt.window.width
      }
    }
	#xlimit: rt.time boundary of scantime
    xlimit <-c(which(xfile@scantime>rt.min)[1]-1, which(xfile@scantime>rt.max)[1] )
    if (is.na(xlimit[2]) ) xlimit[2] <- length(xfile@scantime)
	#ylimit: intensity range of the ion
    ylimit <- range(c(raw.ECI.light[[2]][xlimit[1]:xlimit[2]], raw.ECI.heavy[[2]][xlimit[1]:xlimit[2]], raw.ECI.medium[[2]][xlimit[1]:xlimit[2]]))
    ylimit[1] <- 0.0
    ylimit[2] <- ylimit[2]*1.2
    local.xlimit <- xlimit <- c(rt.min,rt.max)/60
    raw.ECI.light.rt <- xfile@scantime[ raw.ECI.light[[1]] ] / 60
	raw.ECI.medium.rt <- xfile@scantime[ raw.ECI.medium[[1]] ] / 60
    raw.ECI.heavy.rt <- xfile@scantime[ raw.ECI.heavy[[1]] ] / 60
    #title of the EIC profile
    tt.main <- paste(tag, "L/H",raw.file, "; Raw Scan:", as.character(raw.scan.num[j]),
                     "; NL:", formatC(ylimit[2], digits=2, format="e"))
	#plot EIC picture
    plot(raw.ECI.light.rt, raw.ECI.light[[2]], type="l", col="red",xlab="Retention Time(min)",
         ylab="intensity", main=tt.main, xlim=xlimit,ylim=ylimit)
    lines(raw.ECI.heavy.rt, raw.ECI.heavy[[2]], col='blue', xlim=xlimit, ylim=ylimit)
	##the vector containing all of the ms1 scan num and rt time of this peptide
    k_ms1_rt_v <- k_ms1_scan.v <- numeric(0)
    k.ms1.int.light.v <- k_ms1_int_heavy_v <- 0
    if ( !is.na(tag.rt) ) {
	#the ms2 scan num list of this key(this peptide)
      all.ms2.scan <- as.integer( all.scan.table[(key==all.scan.table[,"key"]
                                                  &cross.vec[j]==all.scan.table[,"run"]),"scan"] )
	#the heavy or light list of this key(this peptide)
      all.ms2.HL <- all.scan.table[(key==all.scan.table[,"key"]
                                    &cross.vec[j]==all.scan.table[,"run"]),"HL"]
	#this for  loop is to point the intensity of the light or heavy peptide on the EIC picture
      for (k in 1:length(all.ms2.scan)) {
		#ms1 scan num of this ms2 scan
        k_ms1_scan <- which(xfile@acquisitionNum > all.ms2.scan[k])[1]-1
		#if this ms1 scan num does not exit
        if (is.na(k_ms1_scan)) {
          k_ms1_scan <- length(xfile@acquisitionNum)
        }
		#rt time of the ms1 scan
        k_ms1_rt <- xfile@scantime[k_ms1_scan]/60
		#
        if (all.ms2.HL[k] == "light") {
          points(k_ms1_rt, raw.ECI.light[[2]][k_ms1_scan], type='p',cex=0.5, pch=1)
          #k.ms1.int.light.v <- c(k.ms1.int.light.v, raw.ECI.light[[2]][k_ms1_scan])
        } else if (all.ms2.HL[k] == "heavy") {
          points(k_ms1_rt, raw.ECI.heavy[[2]][k_ms1_scan], type='p',cex=0.5, pch=1)
          #k_ms1_int_heavy_v <- c(k_ms1_int_heavy_v, raw.ECI.heavy[[2]][k_ms1_scan])
        } else {
          points(k_ms1_rt, max(raw.ECI.light[[2]][k_ms1_scan],raw.ECI.heavy[[2]][k_ms1_scan]),
                 type='p',cex=0.5, pch=1,col="black")
        }
        k_ms1_rt_v <- c(k_ms1_rt_v,k_ms1_rt)
        k_ms1_scan.v <- c(k_ms1_scan.v,k_ms1_scan)
      }
      ##lines(c(tag.rt,tag.rt),c(0.0, max(raw.ECI.light[[2]],raw.ECI.heavy[[2]])), col="green")
      HL <- all.ms2.HL[k_ms1_scan.v == tag_ms1_scan_num][1]
      if (HL == "light") {
        points(tag.rt, raw.ECI.light[[2]][tag_ms1_scan_num], type='p',pch=8)
      } else if (HL == "heavy") {
        points(tag.rt, raw.ECI.heavy[[2]][tag_ms1_scan_num], type='p',pch=8)
      } else {
        points(tag.rt, max(raw.ECI.light[[2]][tag_ms1_scan_num],raw.ECI.heavy[[2]][tag_ms1_scan_num]), type='p',pch=8)
      }
      ## record MS1 intensity at which MS2 is triggered
      k.ms1.int.light.v <- raw.ECI.light[[2]][tag_ms1_scan_num]
	  k_ms1_int_medium_v <- raw.ECI.medium[[2]][tag_ms1_scan_num]
      k_ms1_int_heavy_v <- raw.ECI.heavy[[2]][tag_ms1_scan_num]
      ## guess ratio of integrated peak area
      local.xlimit <- c(max(scan.time.range[1]/60, tag.rt-local.rt.window),
                        min(scan.time.range[2]/60, tag.rt+local.rt.window))
    }
    ## guess ratio of integrated peak area
	##return light noise and heavy noise inthe first two position, and rt.min and rt.max pair for paired peaks behind 
    ##see commits for findPairchromPeaks in msisotope.R
    peaks <- findPairChromPeaks( raw.ECI.light.rt, raw.ECI.light[[2]], raw.ECI.heavy[[2]],
                                xlimit, local.xlimit, sn )

	correction.factor <- correction.factor.heavy
	mass_shift <- mass_shift.heavy
    noise.light <- peaks[1]
    lines(xlimit,c(noise.light, noise.light), col='red', type='l', lty=2)
    noise.heavy <- peaks[2]
    lines(xlimit,c(noise.heavy, noise.heavy), col='blue', type='l', lty=2)

    #delete the noise.light and heavy
    peaks <- peaks[-c(1,2)]
    #valid paired peak num
    n.peaks <- length(peaks)/2

    best.yes.single <- best.peak.scan.num <- best.mono_check <- best.r2 <- best.npoints <- best.light.int <- best.ratio <- 0.0
    best.mono_check <- -0.1
    best_xlm <- best.light.yes <- best.heavy.yes <- best.low <- best.high <- c(0)
    best_fixed <- F
    n.light.ms2.peak <- n.heavy.ms2.peak <- n.candidate.peaks <- n_ms2_peaks <- 0
    if (n.peaks>0) {
      for (n in 1:n.peaks) {
        low <- peaks[2*n-1]
        high<- peaks[2*n]
        ### when requested, choose peaks with ms2 events only ###
        if (peaks.with.ms2.only | !is.na(tag.rt)) {
          if (length(k_ms1_rt_v>0) & (sum((k_ms1_rt_v>=low & k_ms1_rt_v<=high))<=0)) next
        }
        yes <- which( raw.ECI.light.rt>=low & raw.ECI.light.rt<=high )
        light.yes <- raw.ECI.light[[2]][yes]
        heavy.yes <- raw.ECI.heavy[[2]][yes]

        peak.scan.num <- raw.ECI.light[[1]][yes][which.max(light.yes)]
        if ( mass_shift > 0 ) {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.heavy) )
        } else {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.heavy+5) )
        }
        #checkChargeAndMonoMass takes the predicted.dist and the peak.scan as the input, 
		#output the correlation factor between the predicted distribution and the experimental distribution
        mono_check <- checkChargeAndMonoMass( peak.scan, mono.mass, charge, mz.ppm.cut, predicted.dist)
		peak.scan.num.heavy <- raw.ECI.heavy[[1]][yes][which.max(heavy.yes)]
        peak.scan.heavy <- getScan(xfile, peak.scan.num.heavy, mzrange=c((mono.mass.heavy-2)/charge, mz.heavy+5) )
		mono_check.heavy <- checkChargeAndMonoMass( peak.scan.heavy, mono.mass.heavy, charge, mz.ppm.cut, predicted.dist.heavy[(mass_shift+1):length(predicted.dist.heavy)])
		mono_check <- max(mono_check, mono_check.heavy)
		if (mono_check == mono_check.heavy) {
		  peak.scan <- peak.scan.heavy
		  peak.scan.num <- peak.scan.num.heavy
		}
        ## calculate ratio of integrated peak area
        ## if we want the H/L ratio, we need to calculate the mono_check for heavy and compare it with the mono_check light, finall take the max one as the mono_check
        if (HL.ratios[j]) {
          ratio <- round((sum(heavy.yes)/sum(light.yes))*correction.factor,digits=2)
        } else {
          ratio <- round((sum(light.yes)/sum(heavy.yes))/correction.factor,digits=2)
        }
		#e.g. in L/H, whether L/H > singletong.ratio or L/H < singleton.ratio, let singleton checker handle this case
        if( singleton.ratio > 0  & ( ratio > singleton.ratio | ratio < 1/singleton.ratio )) next ## let singleton checker handle this case
        lines(c(low,low),ylimit/10, col="green")
        lines(c(high,high),ylimit/10, col="green")
        text(mean(c(low,high)),max(light.yes,heavy.yes)*1.2,
             labels=paste(round(ratio,2),round(mono_check,2),sep="/"))
        ## calculate peak co-elution profile using only points above noise line
        ##yes2 <- light.yes > noise.light & heavy.yes > noise.heavy
        ##light.yes <- light.yes[yes2]
        ##heavy.yes <- heavy.yes[yes2]
        if (ratio > ratio.range[2] | ratio < ratio.range[1]) next
        ## peaks not passing envelope score filter
        ##if (mono_check < env.score.cutoff) next ## disable env.score.cutoff check here as it will be handled by cimage_combine exlusively.

        ## peaks are too narrow
        npoints <- length(light.yes)
        if (npoints<minimum.peak.points) {
          next
        }
        ## extra information for better filtering
        #n_ms2_peaks: the ms1 scan num in this rt window and having ms2 detected
        if (length(k_ms1_rt_v>0) & (sum((k_ms1_rt_v>=low & k_ms1_rt_v<=high))>0)) {
          n_ms2_peaks <- n_ms2_peaks + 1
        }
        x.lm <- lsfit( x=heavy.yes, y=light.yes,intercept=F )
        r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
        #if statement to determine wheter this peak is valid
        if (r2>r2.cutoff) {
          n.candidate.peaks <- n.candidate.peaks + 1
        }
        #if this tag.rt is in this peak rt window, this is best fixed peak
        if ( !is.na(tag.rt) & tag.rt>=low & tag.rt<=high) {
          best_fixed <- T
        } else {
          best_fixed <- F
        }
        if ( best_fixed | (best.mono_check < 0.95 & mono_check >= best.mono_check) | ## better envelope score
            ( best.mono_check >=0.95 & mono_check >=0.95 & max(light.yes,heavy.yes)>max(best.light.yes) ) | ## envelope score equally better, choose a higher peak
            ( mono_check < env.score.cutoff & mono_check == best.mono_check & max(light.yes,heavy.yes) > max(best.light.yes) ) ## envelope score equally worse, choose a higher peak
            ) {
          best.mono_check <- mono_check
          best.npoints <- npoints
          best.r2 <- r2
          best.ratio <- ratio
          best.light.int <- sum(light.yes)
          best_xlm <- round(as.numeric(ls.print(x.lm,print.it=F)$coef.table[[1]][,"Estimate"]),digits=2)+0.01
          best.low <- low
          best.high <- high
          best.light.yes <- light.yes
          best.heavy.yes <- heavy.yes
          best.peak.scan.num <- peak.scan.num
        }
        #if the best peak has already been fixed, break is for loop
        if (best_fixed) break
      }
    }


    if (!best_fixed & !is.na(tag.rt) & (singleton.ratio>0) & HL != "medium") { # if no MS2 within a peak pair, try to identify a singleton peak with MS2
      singleton_ms2_match <- T # whether a singleton peak has a matching MS2, e.g., light for light or heavy for heavy
	  #using HL.ratios to control the HL singleton or LH singleton
	    #singleton <- "L"
        if (HL == "light") { 
			singleton <- "L"  # for singleton heavy, if ms2 is from light, skip
			#mono.single <- mono.mass.heavy
			#raw.ECI.rt.single <- raw.ECI.heavy.rt
			#raw.ECI.single <- raw.ECI.heavy
			#calculate the ratio of the maximum heavy intensity and the light intensity ratio
			#k.ms1.int.ratio <- max(k_ms1_int_heavy_v)/max(c(0.01,k.ms1.int.light.v))
			#if the H/L ratio is small, calculate the L/H ratio to find the L/H singleton
			#if (k.ms1.int.ratio < 3) {
			#	singleton <- "L"
				mono.single <- mono.mass
				raw.ECI.rt.single <- raw.ECI.light.rt
				raw.ECI.single <- raw.ECI.light
				raw.ECI.single.other <- raw.ECI.heavy
				k.ms1.int.ratio <- max(k.ms1.int.light.v)/max(c(0.01,k_ms1_int_heavy_v))
		}
		if (HL == "heavy") {
			singleton <- "H" 
			#mono.single <- mono.mass
			#raw.ECI.rt.single <- raw.ECI.light.rt
			#raw.ECI.single <- raw.ECI.light
			#raw.ECI.single.other <- raw.ECI.heavy
			#k.ms1.int.ratio <- max(k.ms1.int.light.v)/max(c(0.01,k_ms1_int_heavy_v))
			#if (k.ms1.int.ratio < 3) {
			#	singleton <- "H"
			mono.single <- mono.mass.heavy
			raw.ECI.rt.single <- raw.ECI.heavy.rt
			raw.ECI.single <- raw.ECI.heavy
			raw.ECI.single.other <- raw.ECI.light
			k.ms1.int.ratio <- max(k_ms1_int_heavy_v)/max(c(0.01,k.ms1.int.light.v))
		}
      
      n.single.peaks <- 0
	  #if (singleton_ms2_match) {# find singleton peaks only when a matching MS2 exists
		single.peaks <- findSingleChromPeaks(raw.ECI.rt.single, raw.ECI.single[[2]],xlimit, local.xlimit, sn )
		single.peaks <- single.peaks[-1]
		n.single.peaks <- length(single.peaks)/2
	  #}
	  singleton.fixed <- F
	  n.singleton.peaks <- numeric(0)

      if (n.single.peaks > 0) {
        for (ns in 1:n.single.peaks) {
          low.single <- single.peaks[2*ns-1]
          high.single <- single.peaks[2*ns]
          k_ms1_rt_v.tmp <- (k_ms1_rt_v >=low.single & k_ms1_rt_v <= high.single)
          if (length(k_ms1_rt_v>0) & (sum(k_ms1_rt_v.tmp)<=0)) next ## skip a singleton peak without MS2
          if (k.ms1.int.ratio < 3) next ## ratio is too small
          yes.single <- which( raw.ECI.rt.single>=low.single & raw.ECI.rt.single<=high.single )
          int.yes.single <- raw.ECI.single[[2]][yes.single]
          int.yes.single.other <- raw.ECI.single.other[[2]][yes.single]
          peak.scan.num <- raw.ECI.single[[1]][yes.single][which.max(int.yes.single)]
		  #for a scan defined by peak.scan.num, get all the mz and intensity between specific mz range
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.single-2)/charge,
                                                       (mono.single+10)/charge) )
          mono_check.single <- checkChargeAndMonoMass( peak.scan, mono.single, charge, mz.ppm.cut,
                                                      predicted.dist)
          lines(c(low.single,low.single),ylimit/2,col="green")
          lines(c(high.single,high.single),ylimit/2,col="green")
		  # if we get the light or heavy singleton we want, calculate it as normal ratio
		  if ((!HL.ratios[j] & singleton == "H") | (HL.ratios[j] & singleton == "L")) {
			real.singleton.ratio <- round(max(1/singleton.ratio, sum(int.yes.single.other)/sum(int.yes.single)), 2)
		  }
		  # if we get the opposite singleton (e.g. heavy singleton for light ratio), still calculate it as the normal 
		  if ((!HL.ratios[j] & singleton == "L") | (HL.ratios[j] & singleton == "H")) {
		    real.singleton.ratio <- round(min(singleton.ratio, sum(int.yes.single)/max(1,sum(int.yes.single.other))), 2)
		  }
          text(mean(c(low.single,high.single)),max(int.yes.single)*1.2, labels=paste(round(real.singleton.ratio,2),round(mono_check.single,2),sep="/"))
          #did not pass env score filter
          #if (mono_check.single < env.score.cutoff) next
          npoints.single <- length(yes.single)
          ## peak is too narrow with very few time points
          if (npoints.single<minimum.peak.points) next
          #i_ratios[j] <- singleton.ratio
          #r2.v[j] <- 1.0
		  if ( !is.na(tag.rt) & tag.rt>=low.single & tag.rt<=high.single) {
			singleton.fixed <- T
		  } else {
			singleton.fixed <- F
		  }
		  if ( singleton.fixed | (best.mono_check < 0.95 & mono_check.single >= best.mono_check) | ## better envelope score
            ( best.mono_check >=0.95 & mono_check.single >=0.95 & max(int.yes.single) > max(best.yes.single))  | ## envelope score equally better, choose a higher peak
            ( mono_check.single < env.score.cutoff & mono_check.single == best.mono_check & max(int.yes.single) > max(best.yes.single) ) ## envelope score equally worse, choose a higher peak
            ) {
			best.mono_check <- mono_check.single
			best.npoints <- npoints.single
			best.r2 <- 1.0
			best.ratio <- real.singleton.ratio ##min(singleton.ratio, sum(int.yes.single)/max(1,sum(int.yes.single.other)))
			best.light.int <- sum(int.yes.single)
			best_xlm <- 0
			best.low <- low.single
			best.high <- high.single
			best.yes.single <- int.yes.single
			best.light.yes <- 0
			best.heavy.yes <- 0
			best.peak.scan.num <- peak.scan.num
			}


          #plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
          n.singleton.peaks <- c(n.singleton.peaks,ns)
          n_ms2_peaks <- n_ms2_peaks + 1
          n.candidate.peaks <- n.candidate.peaks + 1
		  
          if (singleton.fixed) break
        }
      }
      #   if (length(n.singleton.peaks) == 0) {
      #     plot(0,0,xlab="",ylab="",main=paste("R2 value: 0.00") )
      #     plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
      #   }
    }
	#either paired ratio or singleton ratio exist
    if ( best.r2 != 0 ) {
	# 
      if (best.mono_check >= env.score.cutoff) {
        i_ratios_LH[j] <- best.ratio
        light.int.v[j] <- best.light.int
        l_ratios_LH[j] <- best_xlm
        r2.v_LH[j] <- best.r2
        lines(c(best.low,best.low),ylimit, col="green")
        lines(c(best.high,best.high),ylimit, col="green")
      }
	  #for singleton plot
      if (best_xlm == 0) {
        plot(0,0,xlab="",ylab="",main=paste("Singleton Peak",singleton," Found ! Np = ",npoints.single,sep=""),col.main="red" )
      } else {
	  #for pair ratio plot
        plot(best.heavy.yes,best.light.yes,
             xlab="intensity.heavy", ylab="intensity.light",
             main=paste("X=",format(best_xlm,digits=4),"; R2=",format(best.r2,digits=3),
               "; Np=", best.npoints, sep=""),
             xlim=c(0, max(best.light.yes,best.heavy.yes)),
             ylim=c(0, max(best.light.yes,best.heavy.yes)))
      }
      abline(0,best_xlm)
      abline(0,1,col="grey")
      ## plot raw spectrum
      ##predicted.dist <- predicted.dist[1:20]
      ## upper limit: heavy + 20units ##
      cc <- seq(1,max(which(predicted.dist.heavy>0.01)))
	  #get the merged predicted distribution for heavy and light
      if (HL.ratios[j]) {
        predicted.dist.merge <- (1/best.ratio)*predicted.dist[cc] + predicted.dist.heavy[cc]
      } else {
        predicted.dist.merge <- (best.ratio)*predicted.dist[cc] + predicted.dist.heavy[cc]
      }

      mz.unit <- isotope.mass.unit/charge
      ##predicted.mz <- mono.mz + mz.unit*(seq(1,mass_shift)-1)
      light.index <- which(predicted.dist>0.01)-1
      light.index <- light.index[which(light.index<=mass_shift)]
      predicted.mz <- mono.mz + mz.unit*light.index
      predicted.dist.local <- predicted.dist.merge[light.index+1]
      #predicted.mz.heavy <- mono.mz.heavy + mz.unit*(seq(1, length(predicted.dist.merge)-mass_shift)-1)
      mz.unit.N15 <- isotope.mass.unit.N15/charge
      heavy.index <- which(predicted.dist.heavy>0.01)
      predicted.dist.heavy.local <- predicted.dist.merge[heavy.index]
      heavy.adjustments <- heavy.index <- heavy.index-mass_shift-1
      heavy.adjustments[which(heavy.index<0)] <- mz.unit.N15
      heavy.adjustments[which(heavy.index>=0)] <- mz.unit
      predicted.mz.heavy <- mono.mz.heavy + heavy.adjustments*heavy.index
	# for the predicted distribution ,get those beyond 0.01 and the corresbonding predicted mz
      predicted.mz <- c(predicted.mz, predicted.mz.heavy)

      predicted.dist.merge <- c(predicted.dist.local,predicted.dist.heavy.local)
      n.max <- which.max(predicted.dist.merge)
      predicted.dist.merge <- predicted.dist.merge/predicted.dist.merge[n.max]
		#get the observed intensit according to the predicted mz
      mz.max <- predicted.mz[n.max]
      mass.range <- c(mono.mz-2*mz.unit, mz.heavy+8*mz.unit)
      scan.data <- getScan(xfile, best.peak.scan.num, mzrange=mass.range)
      scan.mz <- scan.data[,"mz"]
      scan.int <- scan.data[,"intensity"]
      observed.int <- predicted.mz
      for ( k in 1:length(observed.int) ) {
        this.mz <- predicted.mz[k]
        mz.diff <- abs(scan.mz-this.mz)/this.mz
        if (min(mz.diff) <= mz.ppm.cut ) {
          observed.int[k] <- scan.int[which.min(mz.diff)]
        } else {
          observed.int[k] <- 0.0
        }
      }
      int.max <- observed.int[n.max]
      if ( max(int.max) > 0.0 ) {
        observed.int <- observed.int / int.max
        scan.int <- scan.int / int.max
      } else {
        scan.int <- scan.int/max(scan.int)
      }
      ylimit2 <- c(0,1.1)
	  ##plot the predicted distribution versus observed distribution picture(right botoom)
	  #plot the background ion intensity, scan.int is for the background
      plot(scan.mz, scan.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="gray")
      par(new=T)
	  #plot the observed distribution
      plot(predicted.mz, observed.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="black")
      if ( mass_shift >0 ){
        light.n <- seq(1,length(light.index))##seq(1,(mass_shift))
        heavy.n <- seq(length(light.index)+1, length(predicted.mz))##seq((mass_shift+1),2*mass_shift)
      } else {
        light.n <- heavy.n <- seq(1,3)
      }
	  #plot the predicted distribution
      par(new=T)
      plot( predicted.mz[light.n], predicted.dist.merge[light.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)
      par(new=T)
      plot( predicted.mz[heavy.n], predicted.dist.merge[heavy.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)

      points(predicted.mz[light.n],rep(0,length(light.n)), pch=23,col="red",bg="white")
      points(predicted.mz[heavy.n],rep(0,length(heavy.n)), pch=24,col="blue",bg="white")
      points(mz.light,0, pch=24,col="red",bg="red")
      points(mz.heavy,0, pch=24,col="blue",bg="blue")
      title( paste("Scan # ", xfile@acquisitionNum[best.peak.scan.num], " @ ",
                   round(xfile@scantime[best.peak.scan.num]/60,1)," min; NL:",
                   formatC(int.max, digits=2,format="e"), sep = ""))
    } else {
      plot(0,0,xlab="",ylab="",main=paste("No quantified peak(s)") )
      plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
      #}
    }
	#cat( paste(i,"success", sep = ":"))
#--------------------------------------light/medium--------------------------------------------------------
    tt.main <- paste(tag, "L/M",raw.file, "; Raw Scan:", as.character(raw.scan.num[j]),
                     "; NL:", formatC(ylimit[2], digits=2, format="e"))
	plot(raw.ECI.light.rt, raw.ECI.light[[2]], type="l", col="red",xlab="Retention Time(min)",
         ylab="intensity", main=tt.main, xlim=xlimit,ylim=ylimit)
	lines(raw.ECI.medium.rt, raw.ECI.medium[[2]], col='green', xlim=xlimit, ylim=ylimit)
	##the vector containing all of the ms1 scan num and rt time of this peptide
    if ( !is.na(tag.rt) ) {
      for (k in 1:length(all.ms2.scan)) {
		#if this ms1 scan num does not exit
        if (all.ms2.HL[k] == "light") {
          points(k_ms1_rt, raw.ECI.light[[2]][k_ms1_scan], type='p',cex=0.5, pch=1)
          #k.ms1.int.light.v <- c(k.ms1.int.light.v, raw.ECI.light[[2]][k_ms1_scan])
        } else if (all.ms2.HL[k] == "medium") {
		  points(k_ms1_rt, raw.ECI.medium[[2]][k_ms1_scan], type='p',cex=0.5, pch=1)
		} else {
          points(k_ms1_rt, max(raw.ECI.light[[2]][k_ms1_scan],raw.ECI.medium[[2]][k_ms1_scan]),
                 type='p',cex=0.5, pch=1,col="black")
        }
      }
      ##lines(c(tag.rt,tag.rt),c(0.0, max(raw.ECI.light[[2]],raw.ECI.heavy[[2]])), col="green")
      if (HL == "light") {
        points(tag.rt, raw.ECI.light[[2]][tag_ms1_scan_num], type='p',pch=8)
      } else if (HL == "medium") {
        points(tag.rt, raw.ECI.medium[[2]][tag_ms1_scan_num], type='p',pch=8)
	  } else {
        points(tag.rt, max(raw.ECI.light[[2]][tag_ms1_scan_num],raw.ECI.medium[[2]][tag_ms1_scan_num]), type='p',pch=8)
      }
    }
	peaks <- findPairChromPeaks( raw.ECI.light.rt, raw.ECI.light[[2]], raw.ECI.medium[[2]],
                                xlimit, local.xlimit, sn )
	noise.medium <- peaks[2]
    lines(xlimit,c(noise.medium, noise.medium), col='green', type='l', lty=2)
	noise.light <- peaks[1]
    lines(xlimit,c(noise.light, noise.light), col='red', type='l', lty=2)
	correction.factor <- correction.factor.medium
	mass_shift <- mass_shift.medium
    #delete the noise.light and medium
    peaks <- peaks[-c(1,2)]
    #valid paired peak num
    n.peaks <- length(peaks)/2

    best.yes.single <- best.peak.scan.num <- best.mono_check <- best.r2 <- best.npoints <- best.light.int <- best.ratio <- 0.0
    best.mono_check <- -0.1
    best_xlm <- best.light.yes <- best.medium.yes <- best.low <- best.high <- c(0)
    best_fixed <- F
    n.light.ms2.peak <- n.medium.ms2.peak <- n.candidate.peaks <- n_ms2_peaks <- 0
    if (n.peaks>0) {
      for (n in 1:n.peaks) {
        low <- peaks[2*n-1]
        high<- peaks[2*n]
        ### when requested, choose peaks with ms2 events only ###
        if (peaks.with.ms2.only | !is.na(tag.rt)) {
          if (length(k_ms1_rt_v>0) & (sum((k_ms1_rt_v>=low & k_ms1_rt_v<=high))<=0)) next
        }
        yes <- which( raw.ECI.light.rt>=low & raw.ECI.light.rt<=high )
        light.yes <- raw.ECI.light[[2]][yes]
        medium.yes <- raw.ECI.medium[[2]][yes]

        peak.scan.num <- raw.ECI.light[[1]][yes][which.max(light.yes)]
        if ( mass_shift > 0 ) {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.medium) )
        } else {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.medium+5) )
        }
        #checkChargeAndMonoMass takes the predicted.dist and the peak.scan as the input, 
		#output the correlation factor between the predicted distribution and the experimental distribution
        mono_check <- checkChargeAndMonoMass( peak.scan, mono.mass, charge, mz.ppm.cut, predicted.dist)
		peak.scan.num.medium <- raw.ECI.medium[[1]][yes][which.max(medium.yes)]
        peak.scan.medium <- getScan(xfile, peak.scan.num.medium, mzrange=c((mono.mass.medium-2)/charge, mz.medium+5) )
		mono_check.medium <- checkChargeAndMonoMass( peak.scan.medium, mono.mass.medium, charge, mz.ppm.cut, predicted.dist.medium[(mass_shift+1):length(predicted.dist.medium)])
		mono_check <- max(mono_check, mono_check.medium)
		if (mono_check == mono_check.medium) {
		  peak.scan <- peak.scan.medium
		  peak.scan.num <- peak.scan.num.medium
		}
        ## calculate ratio of integrated peak area
        ## if we want the H/L ratio, we need to calculate the mono_check for medium and compare it with the mono_check light, finall take the max one as the mono_check
        if (HL.ratios[j]) {
          ratio <- round((sum(medium.yes)/sum(light.yes))*correction.factor,digits=2)
        } else {
          ratio <- round((sum(light.yes)/sum(medium.yes))/correction.factor,digits=2)
        }
		#e.g. in L/H, whether L/H > singletong.ratio or L/H < singleton.ratio, let singleton checker handle this case
        if( singleton.ratio > 0  & ( ratio > singleton.ratio | ratio < 1/singleton.ratio )) next ## let singleton checker handle this case
        lines(c(low,low),ylimit/10, col="green")
        lines(c(high,high),ylimit/10, col="green")
        text(mean(c(low,high)),max(light.yes,medium.yes)*1.2,
             labels=paste(round(ratio,2),round(mono_check,2),sep="/"))
        ## calculate peak co-elution profile using only points above noise line
        ##yes2 <- light.yes > noise.light & medium.yes > noise.medium
        ##light.yes <- light.yes[yes2]
        ##medium.yes <- medium.yes[yes2]
        if (ratio > ratio.range[2] | ratio < ratio.range[1]) next
        ## peaks not passing envelope score filter
        ##if (mono_check < env.score.cutoff) next ## disable env.score.cutoff check here as it will be handled by cimage_combine exlusively.

        ## peaks are too narrow
        npoints <- length(light.yes)
        if (npoints<minimum.peak.points) {
          next
        }
        ## extra information for better filtering
        #n_ms2_peaks: the ms1 scan num in this rt window and having ms2 detected
        if (length(k_ms1_rt_v>0) & (sum((k_ms1_rt_v>=low & k_ms1_rt_v<=high))>0)) {
          n_ms2_peaks <- n_ms2_peaks + 1
        }
        x.lm <- lsfit( x=medium.yes, y=light.yes,intercept=F )
        r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
        #if statement to determine wheter this peak is valid
        if (r2>r2.cutoff) {
          n.candidate.peaks <- n.candidate.peaks + 1
        }
        #if this tag.rt is in this peak rt window, this is best fixed peak
        if ( !is.na(tag.rt) & tag.rt>=low & tag.rt<=high) {
          best_fixed <- T
        } else {
          best_fixed <- F
        }
        if ( best_fixed | (best.mono_check < 0.95 & mono_check >= best.mono_check) | ## better envelope score
            ( best.mono_check >=0.95 & mono_check >=0.95 & max(light.yes,medium.yes)>max(best.light.yes) ) | ## envelope score equally better, choose a higher peak
            ( mono_check < env.score.cutoff & mono_check == best.mono_check & max(light.yes,medium.yes) > max(best.light.yes) ) ## envelope score equally worse, choose a higher peak
            ) {
          best.mono_check <- mono_check
          best.npoints <- npoints
          best.r2 <- r2
          best.ratio <- ratio
          best.light.int <- sum(light.yes)
          best_xlm <- round(as.numeric(ls.print(x.lm,print.it=F)$coef.table[[1]][,"Estimate"]),digits=2)+0.01
          best.low <- low
          best.high <- high
          best.light.yes <- light.yes
          best.medium.yes <- medium.yes
          best.peak.scan.num <- peak.scan.num
        }
        #if the best peak has already been fixed, break is for loop
        if (best_fixed) break
      }
    }


    if (!best_fixed & !is.na(tag.rt) & (singleton.ratio>0) & HL != "heavy") { # if no MS2 within a peak pair, try to identify a singleton peak with MS2
      singleton_ms2_match <- T # whether a singleton peak has a matching MS2, e.g., light for light or medium for medium
	  #using HL.ratios to control the HL singleton or LH singleton
	    #singleton <- "L"
        if (HL == "light") { 
			singleton <- "L"  # for singleton medium, if ms2 is from light, skip
			#mono.single <- mono.mass.medium
			#raw.ECI.rt.single <- raw.ECI.medium.rt
			#raw.ECI.single <- raw.ECI.medium
			#calculate the ratio of the maximum medium intensity and the light intensity ratio
			#k.ms1.int.ratio <- max(k_ms1_int_medium_v)/max(c(0.01,k.ms1.int.light.v))
			#if the H/L ratio is small, calculate the L/H ratio to find the L/H singleton
			#if (k.ms1.int.ratio < 3) {
			#	singleton <- "L"
				mono.single <- mono.mass
				raw.ECI.rt.single <- raw.ECI.light.rt
				raw.ECI.single <- raw.ECI.light
				raw.ECI.single.other <- raw.ECI.medium
				k.ms1.int.ratio <- max(k.ms1.int.light.v)/max(c(0.01,k_ms1_int_medium_v))
		}
		if (HL == "medium") {
			singleton <- "H" 
			#mono.single <- mono.mass
			#raw.ECI.rt.single <- raw.ECI.light.rt
			#raw.ECI.single <- raw.ECI.light
			#raw.ECI.single.other <- raw.ECI.medium
			#k.ms1.int.ratio <- max(k.ms1.int.light.v)/max(c(0.01,k_ms1_int_medium_v))
			#if (k.ms1.int.ratio < 3) {
			#	singleton <- "H"
			mono.single <- mono.mass.medium
			raw.ECI.rt.single <- raw.ECI.medium.rt
			raw.ECI.single <- raw.ECI.medium
			raw.ECI.single.other <- raw.ECI.light
			k.ms1.int.ratio <- max(k_ms1_int_medium_v)/max(c(0.01,k.ms1.int.light.v))
		}
      
      n.single.peaks <- 0
	  #if (singleton_ms2_match) {# find singleton peaks only when a matching MS2 exists
		single.peaks <- findSingleChromPeaks(raw.ECI.rt.single, raw.ECI.single[[2]],xlimit, local.xlimit, sn )
		single.peaks <- single.peaks[-1]
		n.single.peaks <- length(single.peaks)/2
	  #}
	  singleton.fixed <- F
	  n.singleton.peaks <- numeric(0)

      if (n.single.peaks > 0) {
        for (ns in 1:n.single.peaks) {
          low.single <- single.peaks[2*ns-1]
          high.single <- single.peaks[2*ns]
          k_ms1_rt_v.tmp <- (k_ms1_rt_v >=low.single & k_ms1_rt_v <= high.single)
          if (length(k_ms1_rt_v>0) & (sum(k_ms1_rt_v.tmp)<=0)) next ## skip a singleton peak without MS2
          if (k.ms1.int.ratio < 3) next ## ratio is too small
          yes.single <- which( raw.ECI.rt.single>=low.single & raw.ECI.rt.single<=high.single )
          int.yes.single <- raw.ECI.single[[2]][yes.single]
          int.yes.single.other <- raw.ECI.single.other[[2]][yes.single]
          peak.scan.num <- raw.ECI.single[[1]][yes.single][which.max(int.yes.single)]
		  #for a scan defined by peak.scan.num, get all the mz and intensity between specific mz range
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.single-2)/charge,
                                                       (mono.single+10)/charge) )
          mono_check.single <- checkChargeAndMonoMass( peak.scan, mono.single, charge, mz.ppm.cut,
                                                      predicted.dist)
          lines(c(low.single,low.single),ylimit/2,col="green")
          lines(c(high.single,high.single),ylimit/2,col="green")
		  # if we get the light or medium singleton we want, calculate it as normal ratio
		  if ((!HL.ratios[j] & singleton == "H") | (HL.ratios[j] & singleton == "L")) {
			real.singleton.ratio <- round(max(1/singleton.ratio, sum(int.yes.single.other)/sum(int.yes.single)), 2)
		  }
		  # if we get the opposite singleton (e.g. medium singleton for light ratio), still calculate it as the normal 
		  if ((!HL.ratios[j] & singleton == "L") | (HL.ratios[j] & singleton == "H")) {
		    real.singleton.ratio <- round(min(singleton.ratio, sum(int.yes.single)/max(1,sum(int.yes.single.other))), 2)
		  }
          text(mean(c(low.single,high.single)),max(int.yes.single)*1.2, labels=paste(round(real.singleton.ratio,2),round(mono_check.single,2),sep="/"))
          #did not pass env score filter
          #if (mono_check.single < env.score.cutoff) next
          npoints.single <- length(yes.single)
          ## peak is too narrow with very few time points
          if (npoints.single<minimum.peak.points) next
          #i_ratios[j] <- singleton.ratio
          #r2.v[j] <- 1.0
		  if ( !is.na(tag.rt) & tag.rt>=low.single & tag.rt<=high.single) {
			singleton.fixed <- T
		  } else {
			singleton.fixed <- F
		  }
		  if ( singleton.fixed | (best.mono_check < 0.95 & mono_check.single >= best.mono_check) | ## better envelope score
            ( best.mono_check >=0.95 & mono_check.single >=0.95 & max(int.yes.single) > max(best.yes.single))  | ## envelope score equally better, choose a higher peak
            ( mono_check.single < env.score.cutoff & mono_check.single == best.mono_check & max(int.yes.single) > max(best.yes.single) ) ## envelope score equally worse, choose a higher peak
            ) {
			best.mono_check <- mono_check.single
			best.npoints <- npoints.single
			best.r2 <- 1.0
			best.ratio <- real.singleton.ratio ##min(singleton.ratio, sum(int.yes.single)/max(1,sum(int.yes.single.other)))
			best.light.int <- sum(int.yes.single)
			best_xlm <- 0
			best.low <- low.single
			best.high <- high.single
			best.yes.single <- int.yes.single
			best.light.yes <- 0
			best.medium.yes <- 0
			best.peak.scan.num <- peak.scan.num
			}


          #plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
          n.singleton.peaks <- c(n.singleton.peaks,ns)
          n_ms2_peaks <- n_ms2_peaks + 1
          n.candidate.peaks <- n.candidate.peaks + 1
		  
          if (singleton.fixed) break
        }
      }
      #   if (length(n.singleton.peaks) == 0) {
      #     plot(0,0,xlab="",ylab="",main=paste("R2 value: 0.00") )
      #     plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
      #   }
    }
	#either paired ratio or singleton ratio exist
    if ( best.r2 != 0 ) {
	# 
      if (best.mono_check >= env.score.cutoff) {
        i_ratios_LM[j] <- best.ratio
        light.int.v[j] <- best.light.int
        l_ratios_LM[j] <- best_xlm
        r2.v_LM[j] <- best.r2
        lines(c(best.low,best.low),ylimit, col="green")
        lines(c(best.high,best.high),ylimit, col="green")
      }
	  #for singleton plot
      if (best_xlm == 0) {
        plot(0,0,xlab="",ylab="",main=paste("Singleton Peak",singleton," Found ! Np = ",npoints.single,sep=""),col.main="red" )
      } else {
	  #for pair ratio plot
        plot(best.medium.yes,best.light.yes,
             xlab="intensity.medium", ylab="intensity.light",
             main=paste("X=",format(best_xlm,digits=4),"; R2=",format(best.r2,digits=3),
               "; Np=", best.npoints, sep=""),
             xlim=c(0, max(best.light.yes,best.medium.yes)),
             ylim=c(0, max(best.light.yes,best.medium.yes)))
      }
      abline(0,best_xlm)
      abline(0,1,col="grey")
      ## plot raw spectrum
      ##predicted.dist <- predicted.dist[1:20]
      ## upper limit: medium + 20units ##
      cc <- seq(1,max(which(predicted.dist.medium>0.01)))
	  #get the merged predicted distribution for medium and light
      if (HL.ratios[j]) {
        predicted.dist.merge <- (1/best.ratio)*predicted.dist[cc] + predicted.dist.medium[cc]
      } else {
        predicted.dist.merge <- (best.ratio)*predicted.dist[cc] + predicted.dist.medium[cc]
      }

      mz.unit <- isotope.mass.unit/charge
      ##predicted.mz <- mono.mz + mz.unit*(seq(1,mass_shift)-1)
      light.index <- which(predicted.dist>0.01)-1
      light.index <- light.index[which(light.index<=mass_shift)]
      predicted.mz <- mono.mz + mz.unit*light.index
      predicted.dist.local <- predicted.dist.merge[light.index+1]
      #predicted.mz.medium <- mono.mz.medium + mz.unit*(seq(1, length(predicted.dist.merge)-mass_shift)-1)
      mz.unit.N15 <- isotope.mass.unit.N15/charge
      medium.index <- which(predicted.dist.medium>0.01)
      predicted.dist.medium.local <- predicted.dist.merge[medium.index]
      medium.adjustments <- medium.index <- medium.index-mass_shift-1
      medium.adjustments[which(medium.index<0)] <- mz.unit.N15
      medium.adjustments[which(medium.index>=0)] <- mz.unit
      predicted.mz.medium <- mono.mz.medium + medium.adjustments*medium.index
	# for the predicted distribution ,get those beyond 0.01 and the corresbonding predicted mz
      predicted.mz <- c(predicted.mz, predicted.mz.medium)

      predicted.dist.merge <- c(predicted.dist.local,predicted.dist.medium.local)
      n.max <- which.max(predicted.dist.merge)
      predicted.dist.merge <- predicted.dist.merge/predicted.dist.merge[n.max]
		#get the observed intensit according to the predicted mz
      mz.max <- predicted.mz[n.max]
      mass.range <- c(mono.mz-2*mz.unit, mz.medium+8*mz.unit)
      scan.data <- getScan(xfile, best.peak.scan.num, mzrange=mass.range)
      scan.mz <- scan.data[,"mz"]
      scan.int <- scan.data[,"intensity"]
      observed.int <- predicted.mz
      for ( k in 1:length(observed.int) ) {
        this.mz <- predicted.mz[k]
        mz.diff <- abs(scan.mz-this.mz)/this.mz
        if (min(mz.diff) <= mz.ppm.cut ) {
          observed.int[k] <- scan.int[which.min(mz.diff)]
        } else {
          observed.int[k] <- 0.0
        }
      }
      int.max <- observed.int[n.max]
      if ( max(int.max) > 0.0 ) {
        observed.int <- observed.int / int.max
        scan.int <- scan.int / int.max
      } else {
        scan.int <- scan.int/max(scan.int)
      }
      ylimit2 <- c(0,1.1)
	  ##plot the predicted distribution versus observed distribution picture(right botoom)
	  #plot the background ion intensity, scan.int is for the background
      plot(scan.mz, scan.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="gray")
      par(new=T)
	  #plot the observed distribution
      plot(predicted.mz, observed.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="black")
      if ( mass_shift >0 ){
        light.n <- seq(1,length(light.index))##seq(1,(mass_shift))
        medium.n <- seq(length(light.index)+1, length(predicted.mz))##seq((mass_shift+1),2*mass_shift)
      } else {
        light.n <- medium.n <- seq(1,3)
      }
	  #plot the predicted distribution
      par(new=T)
      plot( predicted.mz[light.n], predicted.dist.merge[light.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)
      par(new=T)
      plot( predicted.mz[medium.n], predicted.dist.merge[medium.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)

      points(predicted.mz[light.n],rep(0,length(light.n)), pch=23,col="red",bg="white")
      points(predicted.mz[medium.n],rep(0,length(medium.n)), pch=24,col="blue",bg="white")
      points(mz.light,0, pch=24,col="red",bg="red")
      points(mz.medium,0, pch=24,col="blue",bg="blue")
      title( paste("Scan # ", xfile@acquisitionNum[best.peak.scan.num], " @ ",
                   round(xfile@scantime[best.peak.scan.num]/60,1)," min; NL:",
                   formatC(int.max, digits=2,format="e"), sep = ""))
    } else {
      plot(0,0,xlab="",ylab="",main=paste("No quantified peak(s)") )
      plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
      #}
    }
	#cat( paste(i,"success", sep = ":"))
#--------------------------------------------------------light/medium--------------------------------------------------------
#--------------------------------------------------------heavy/medium--------------------------------------------------------
    tt.main <- paste(tag, "H/M",raw.file, "; Raw Scan:", as.character(raw.scan.num[j]),
                     "; NL:", formatC(ylimit[2], digits=2, format="e"))
	plot(raw.ECI.heavy.rt, raw.ECI.heavy[[2]], type="l", col="red",xlab="Retention Time(min)",
         ylab="intensity", main=tt.main, xlim=xlimit,ylim=ylimit)
	lines(raw.ECI.medium.rt, raw.ECI.medium[[2]], col='green', xlim=xlimit, ylim=ylimit)
	##the vector containing all of the ms1 scan num and rt time of this peptide
    if ( !is.na(tag.rt) ) {
      for (k in 1:length(all.ms2.scan)) {
		#if this ms1 scan num does not exit
        if (all.ms2.HL[k] == "heavy") {
          points(k_ms1_rt, raw.ECI.heavy[[2]][k_ms1_scan], type='p',cex=0.5, pch=1)
          #k_ms1_int_heavy_v <- c(k_ms1_int_heavy_v, raw.ECI.heavy[[2]][k_ms1_scan])
        } else if (all.ms2.HL[k] == "medium") {
		  points(k_ms1_rt, raw.ECI.medium[[2]][k_ms1_scan], type='p',cex=0.5, pch=1)
		} else {
          points(k_ms1_rt, max(raw.ECI.heavy[[2]][k_ms1_scan],raw.ECI.medium[[2]][k_ms1_scan]),
                 type='p',cex=0.5, pch=1,col="black")
        }
      }
      ##lines(c(tag.rt,tag.rt),c(0.0, max(raw.ECI.heavy[[2]],raw.ECI.heavy[[2]])), col="green")
      if (HL == "heavy") {
        points(tag.rt, raw.ECI.heavy[[2]][tag_ms1_scan_num], type='p',pch=8)
      } else if (HL == "medium") {
        points(tag.rt, raw.ECI.medium[[2]][tag_ms1_scan_num], type='p',pch=8)
	  } else {
        points(tag.rt, max(raw.ECI.heavy[[2]][tag_ms1_scan_num],raw.ECI.medium[[2]][tag_ms1_scan_num]), type='p',pch=8)
      }
    }
	peaks <- findPairChromPeaks( raw.ECI.heavy.rt, raw.ECI.heavy[[2]], raw.ECI.medium[[2]],
                                xlimit, local.xlimit, sn )
	noise.medium <- peaks[2]
    lines(xlimit,c(noise.medium, noise.medium), col='green', type='l', lty=2)
	noise.heavy <- peaks[1]
    lines(xlimit,c(noise.heavy, noise.heavy), col='red', type='l', lty=2)
	correction.factor <- correction.factor.medium
	mass_shift <- mass_shift.medium
    #delete the noise.heavy and medium
    peaks <- peaks[-c(1,2)]
    #valid paired peak num
    n.peaks <- length(peaks)/2

    best.yes.single <- best.peak.scan.num <- best.mono_check <- best.r2 <- best.npoints <- best.heavy.int <- best.ratio <- 0.0
    best.mono_check <- -0.1
    best_xlm <- best.heavy.yes <- best.medium.yes <- best.low <- best.high <- c(0)
    best_fixed <- F
    n.heavy.ms2.peak <- n.medium.ms2.peak <- n.candidate.peaks <- n_ms2_peaks <- 0
    if (n.peaks>0) {
      for (n in 1:n.peaks) {
        low <- peaks[2*n-1]
        high<- peaks[2*n]
        ### when requested, choose peaks with ms2 events only ###
        if (peaks.with.ms2.only | !is.na(tag.rt)) {
          if (length(k_ms1_rt_v>0) & (sum((k_ms1_rt_v>=low & k_ms1_rt_v<=high))<=0)) next
        }
        yes <- which( raw.ECI.heavy.rt>=low & raw.ECI.heavy.rt<=high )
        heavy.yes <- raw.ECI.heavy[[2]][yes]
        medium.yes <- raw.ECI.medium[[2]][yes]

        peak.scan.num <- raw.ECI.heavy[[1]][yes][which.max(heavy.yes)]
        if ( mass_shift > 0 ) {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.medium) )
        } else {
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.mass-2)/charge, mz.medium+5) )
        }
        #checkChargeAndMonoMass takes the predicted.dist and the peak.scan as the input, 
		#output the correlation factor between the predicted distribution and the experimental distribution
        mono_check <- checkChargeAndMonoMass( peak.scan, mono.mass, charge, mz.ppm.cut, predicted.dist)
		peak.scan.num.medium <- raw.ECI.medium[[1]][yes][which.max(medium.yes)]
        peak.scan.medium <- getScan(xfile, peak.scan.num.medium, mzrange=c((mono.mass.medium-2)/charge, mz.medium+5) )
		mono_check.medium <- checkChargeAndMonoMass( peak.scan.medium, mono.mass.medium, charge, mz.ppm.cut, predicted.dist.medium[(mass_shift+1):length(predicted.dist.medium)])
		mono_check <- max(mono_check, mono_check.medium)
		if (mono_check == mono_check.medium) {
		  peak.scan <- peak.scan.medium
		  peak.scan.num <- peak.scan.num.medium
		}
        ## calculate ratio of integrated peak area
        ## if we want the H/L ratio, we need to calculate the mono_check for medium and compare it with the mono_check heavy, finall take the max one as the mono_check
        if (HL.ratios[j]) {
          ratio <- round((sum(medium.yes)/sum(heavy.yes))*correction.factor,digits=2)
        } else {
          ratio <- round((sum(heavy.yes)/sum(medium.yes))/correction.factor,digits=2)
        }
		#e.g. in L/H, whether L/H > singletong.ratio or L/H < singleton.ratio, let singleton checker handle this case
        if( singleton.ratio > 0  & ( ratio > singleton.ratio | ratio < 1/singleton.ratio )) next ## let singleton checker handle this case
        lines(c(low,low),ylimit/10, col="green")
        lines(c(high,high),ylimit/10, col="green")
        text(mean(c(low,high)),max(heavy.yes,medium.yes)*1.2,
             labels=paste(round(ratio,2),round(mono_check,2),sep="/"))
        ## calculate peak co-elution profile using only points above noise line
        ##yes2 <- heavy.yes > noise.heavy & medium.yes > noise.medium
        ##heavy.yes <- heavy.yes[yes2]
        ##medium.yes <- medium.yes[yes2]
        if (ratio > ratio.range[2] | ratio < ratio.range[1]) next
        ## peaks not passing envelope score filter
        ##if (mono_check < env.score.cutoff) next ## disable env.score.cutoff check here as it will be handled by cimage_combine exlusively.

        ## peaks are too narrow
        npoints <- length(heavy.yes)
        if (npoints<minimum.peak.points) {
          next
        }
        ## extra information for better filtering
        #n_ms2_peaks: the ms1 scan num in this rt window and having ms2 detected
        if (length(k_ms1_rt_v>0) & (sum((k_ms1_rt_v>=low & k_ms1_rt_v<=high))>0)) {
          n_ms2_peaks <- n_ms2_peaks + 1
        }
        x.lm <- lsfit( x=medium.yes, y=heavy.yes,intercept=F )
        r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
        #if statement to determine wheter this peak is valid
        if (r2>r2.cutoff) {
          n.candidate.peaks <- n.candidate.peaks + 1
        }
        #if this tag.rt is in this peak rt window, this is best fixed peak
        if ( !is.na(tag.rt) & tag.rt>=low & tag.rt<=high) {
          best_fixed <- T
        } else {
          best_fixed <- F
        }
        if ( best_fixed | (best.mono_check < 0.95 & mono_check >= best.mono_check) | ## better envelope score
            ( best.mono_check >=0.95 & mono_check >=0.95 & max(heavy.yes,medium.yes)>max(best.heavy.yes) ) | ## envelope score equally better, choose a higher peak
            ( mono_check < env.score.cutoff & mono_check == best.mono_check & max(heavy.yes,medium.yes) > max(best.heavy.yes) ) ## envelope score equally worse, choose a higher peak
            ) {
          best.mono_check <- mono_check
          best.npoints <- npoints
          best.r2 <- r2
          best.ratio <- ratio
          best.heavy.int <- sum(heavy.yes)
          best_xlm <- round(as.numeric(ls.print(x.lm,print.it=F)$coef.table[[1]][,"Estimate"]),digits=2)+0.01
          best.low <- low
          best.high <- high
          best.heavy.yes <- heavy.yes
          best.medium.yes <- medium.yes
          best.peak.scan.num <- peak.scan.num
        }
        #if the best peak has already been fixed, break is for loop
        if (best_fixed) break
      }
    }


    if (!best_fixed & !is.na(tag.rt) & (singleton.ratio>0) & HL != "light") { # if no MS2 within a peak pair, try to identify a singleton peak with MS2
      singleton_ms2_match <- T # whether a singleton peak has a matching MS2, e.g., heavy for heavy or medium for medium
	  #using HL.ratios to control the HL singleton or LH singleton
	    #singleton <- "L"
        if (HL == "heavy") { 
			singleton <- "L"  # for singleton medium, if ms2 is from heavy, skip
			#mono.single <- mono.mass.medium
			#raw.ECI.rt.single <- raw.ECI.medium.rt
			#raw.ECI.single <- raw.ECI.medium
			#calculate the ratio of the maximum medium intensity and the heavy intensity ratio
			#k.ms1.int.ratio <- max(k_ms1_int_medium_v)/max(c(0.01,k_ms1_int_heavy_v))
			#if the H/L ratio is small, calculate the L/H ratio to find the L/H singleton
			#if (k.ms1.int.ratio < 3) {
			#	singleton <- "L"
				mono.single <- mono.mass
				raw.ECI.rt.single <- raw.ECI.heavy.rt
				raw.ECI.single <- raw.ECI.heavy
				raw.ECI.single.other <- raw.ECI.medium
				k.ms1.int.ratio <- max(k_ms1_int_heavy_v)/max(c(0.01,k_ms1_int_medium_v))
		}
		if (HL == "medium") {
			singleton <- "H" 
			#mono.single <- mono.mass
			#raw.ECI.rt.single <- raw.ECI.heavy.rt
			#raw.ECI.single <- raw.ECI.heavy
			#raw.ECI.single.other <- raw.ECI.medium
			#k.ms1.int.ratio <- max(k_ms1_int_heavy_v)/max(c(0.01,k_ms1_int_medium_v))
			#if (k.ms1.int.ratio < 3) {
			#	singleton <- "H"
			mono.single <- mono.mass.medium
			raw.ECI.rt.single <- raw.ECI.medium.rt
			raw.ECI.single <- raw.ECI.medium
			raw.ECI.single.other <- raw.ECI.heavy
			k.ms1.int.ratio <- max(k_ms1_int_medium_v)/max(c(0.01,k_ms1_int_heavy_v))
		}
      
      n.single.peaks <- 0
	  #if (singleton_ms2_match) {# find singleton peaks only when a matching MS2 exists
		single.peaks <- findSingleChromPeaks(raw.ECI.rt.single, raw.ECI.single[[2]],xlimit, local.xlimit, sn )
		single.peaks <- single.peaks[-1]
		n.single.peaks <- length(single.peaks)/2
	  #}
	  singleton.fixed <- F
	  n.singleton.peaks <- numeric(0)

      if (n.single.peaks > 0) {
        for (ns in 1:n.single.peaks) {
          low.single <- single.peaks[2*ns-1]
          high.single <- single.peaks[2*ns]
          k_ms1_rt_v.tmp <- (k_ms1_rt_v >=low.single & k_ms1_rt_v <= high.single)
          if (length(k_ms1_rt_v>0) & (sum(k_ms1_rt_v.tmp)<=0)) next ## skip a singleton peak without MS2
          if (k.ms1.int.ratio < 3) next ## ratio is too small
          yes.single <- which( raw.ECI.rt.single>=low.single & raw.ECI.rt.single<=high.single )
          int.yes.single <- raw.ECI.single[[2]][yes.single]
          int.yes.single.other <- raw.ECI.single.other[[2]][yes.single]
          peak.scan.num <- raw.ECI.single[[1]][yes.single][which.max(int.yes.single)]
		  #for a scan defined by peak.scan.num, get all the mz and intensity between specific mz range
          peak.scan <- getScan(xfile, peak.scan.num, mzrange=c((mono.single-2)/charge,
                                                       (mono.single+10)/charge) )
          mono_check.single <- checkChargeAndMonoMass( peak.scan, mono.single, charge, mz.ppm.cut,
                                                      predicted.dist)
          lines(c(low.single,low.single),ylimit/2,col="green")
          lines(c(high.single,high.single),ylimit/2,col="green")
		  # if we get the heavy or medium singleton we want, calculate it as normal ratio
		  if ((!HL.ratios[j] & singleton == "H") | (HL.ratios[j] & singleton == "L")) {
			real.singleton.ratio <- round(max(1/singleton.ratio, sum(int.yes.single.other)/sum(int.yes.single)), 2)
		  }
		  # if we get the opposite singleton (e.g. medium singleton for heavy ratio), still calculate it as the normal 
		  if ((!HL.ratios[j] & singleton == "L") | (HL.ratios[j] & singleton == "H")) {
		    real.singleton.ratio <- round(min(singleton.ratio, sum(int.yes.single)/max(1,sum(int.yes.single.other))), 2)
		  }
          text(mean(c(low.single,high.single)),max(int.yes.single)*1.2, labels=paste(round(real.singleton.ratio,2),round(mono_check.single,2),sep="/"))
          #did not pass env score filter
          #if (mono_check.single < env.score.cutoff) next
          npoints.single <- length(yes.single)
          ## peak is too narrow with very few time points
          if (npoints.single<minimum.peak.points) next
          #i_ratios[j] <- singleton.ratio
          #r2.v[j] <- 1.0
		  if ( !is.na(tag.rt) & tag.rt>=low.single & tag.rt<=high.single) {
			singleton.fixed <- T
		  } else {
			singleton.fixed <- F
		  }
		  if ( singleton.fixed | (best.mono_check < 0.95 & mono_check.single >= best.mono_check) | ## better envelope score
            ( best.mono_check >=0.95 & mono_check.single >=0.95 & max(int.yes.single) > max(best.yes.single))  | ## envelope score equally better, choose a higher peak
            ( mono_check.single < env.score.cutoff & mono_check.single == best.mono_check & max(int.yes.single) > max(best.yes.single) ) ## envelope score equally worse, choose a higher peak
            ) {
			best.mono_check <- mono_check.single
			best.npoints <- npoints.single
			best.r2 <- 1.0
			best.ratio <- real.singleton.ratio ##min(singleton.ratio, sum(int.yes.single)/max(1,sum(int.yes.single.other)))
			best.heavy.int <- sum(int.yes.single)
			best_xlm <- 0
			best.low <- low.single
			best.high <- high.single
			best.yes.single <- int.yes.single
			best.heavy.yes <- 0
			best.medium.yes <- 0
			best.peak.scan.num <- peak.scan.num
			}


          #plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
          n.singleton.peaks <- c(n.singleton.peaks,ns)
          n_ms2_peaks <- n_ms2_peaks + 1
          n.candidate.peaks <- n.candidate.peaks + 1
		  
          if (singleton.fixed) break
        }
      }
      #   if (length(n.singleton.peaks) == 0) {
      #     plot(0,0,xlab="",ylab="",main=paste("R2 value: 0.00") )
      #     plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
      #   }
    }
	#either paired ratio or singleton ratio exist
    if ( best.r2 != 0 ) {
	# 
      if (best.mono_check >= env.score.cutoff) {
        i_ratios_HM[j] <- best.ratio
        light.int.v[j] <- best.heavy.int
        l_ratios_HM[j] <- best_xlm
        r2.v_HM[j] <- best.r2
        lines(c(best.low,best.low),ylimit, col="green")
        lines(c(best.high,best.high),ylimit, col="green")
      }
	  #for singleton plot
      if (best_xlm == 0) {
        plot(0,0,xlab="",ylab="",main=paste("Singleton Peak",singleton," Found ! Np = ",npoints.single,sep=""),col.main="red" )
      } else {
	  #for pair ratio plot
        plot(best.medium.yes,best.heavy.yes,
             xlab="intensity.medium", ylab="intensity.heavy",
             main=paste("X=",format(best_xlm,digits=4),"; R2=",format(best.r2,digits=3),
               "; Np=", best.npoints, sep=""),
             xlim=c(0, max(best.heavy.yes,best.medium.yes)),
             ylim=c(0, max(best.heavy.yes,best.medium.yes)))
      }
      abline(0,best_xlm)
      abline(0,1,col="grey")
      ## plot raw spectrum
      ##predicted.dist <- predicted.dist[1:20]
      ## upper limit: medium + 20units ##
      cc <- seq(1,max(which(predicted.dist.medium>0.01)))
	  #get the merged predicted distribution for medium and heavy
      if (HL.ratios[j]) {
        predicted.dist.merge <- (1/best.ratio)*predicted.dist[cc] + predicted.dist.medium[cc]
      } else {
        predicted.dist.merge <- (best.ratio)*predicted.dist[cc] + predicted.dist.medium[cc]
      }

      mz.unit <- isotope.mass.unit/charge
      ##predicted.mz <- mono.mz + mz.unit*(seq(1,mass_shift)-1)
      heavy.index <- which(predicted.dist>0.01)-1
      heavy.index <- heavy.index[which(heavy.index<=mass_shift)]
      predicted.mz <- mono.mz + mz.unit*heavy.index
      predicted.dist.local <- predicted.dist.merge[heavy.index+1]
      #predicted.mz.medium <- mono.mz.medium + mz.unit*(seq(1, length(predicted.dist.merge)-mass_shift)-1)
      mz.unit.N15 <- isotope.mass.unit.N15/charge
      medium.index <- which(predicted.dist.medium>0.01)
      predicted.dist.medium.local <- predicted.dist.merge[medium.index]
      medium.adjustments <- medium.index <- medium.index-mass_shift-1
      medium.adjustments[which(medium.index<0)] <- mz.unit.N15
      medium.adjustments[which(medium.index>=0)] <- mz.unit
      predicted.mz.medium <- mono.mz.medium + medium.adjustments*medium.index
	# for the predicted distribution ,get those beyond 0.01 and the corresbonding predicted mz
      predicted.mz <- c(predicted.mz, predicted.mz.medium)

      predicted.dist.merge <- c(predicted.dist.local,predicted.dist.medium.local)
      n.max <- which.max(predicted.dist.merge)
      predicted.dist.merge <- predicted.dist.merge/predicted.dist.merge[n.max]
		#get the observed intensit according to the predicted mz
      mz.max <- predicted.mz[n.max]
      mass.range <- c(mono.mz-2*mz.unit, mz.medium+8*mz.unit)
      scan.data <- getScan(xfile, best.peak.scan.num, mzrange=mass.range)
      scan.mz <- scan.data[,"mz"]
      scan.int <- scan.data[,"intensity"]
      observed.int <- predicted.mz
      for ( k in 1:length(observed.int) ) {
        this.mz <- predicted.mz[k]
        mz.diff <- abs(scan.mz-this.mz)/this.mz
        if (min(mz.diff) <= mz.ppm.cut ) {
          observed.int[k] <- scan.int[which.min(mz.diff)]
        } else {
          observed.int[k] <- 0.0
        }
      }
      int.max <- observed.int[n.max]
      if ( max(int.max) > 0.0 ) {
        observed.int <- observed.int / int.max
        scan.int <- scan.int / int.max
      } else {
        scan.int <- scan.int/max(scan.int)
      }
      ylimit2 <- c(0,1.1)
	  ##plot the predicted distribution versus observed distribution picture(right botoom)
	  #plot the background ion intensity, scan.int is for the background
      plot(scan.mz, scan.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="gray")
      par(new=T)
	  #plot the observed distribution
      plot(predicted.mz, observed.int, type='h', xlab="m/z", ylab="intensity", xlim=mass.range, ylim=ylimit2, col="black")
      if ( mass_shift >0 ){
        heavy.n <- seq(1,length(heavy.index))##seq(1,(mass_shift))
        medium.n <- seq(length(heavy.index)+1, length(predicted.mz))##seq((mass_shift+1),2*mass_shift)
      } else {
        heavy.n <- medium.n <- seq(1,3)
      }
	  #plot the predicted distribution
      par(new=T)
      plot( predicted.mz[heavy.n], predicted.dist.merge[heavy.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)
      par(new=T)
      plot( predicted.mz[medium.n], predicted.dist.merge[medium.n], type='b',xlab="",ylab="",col="green",axes=F,xlim=mass.range,ylim=ylimit2)

      points(predicted.mz[heavy.n],rep(0,length(heavy.n)), pch=23,col="red",bg="white")
      points(predicted.mz[medium.n],rep(0,length(medium.n)), pch=24,col="blue",bg="white")
      points(mz.heavy,0, pch=24,col="red",bg="red")
      points(mz.medium,0, pch=24,col="blue",bg="blue")
      title( paste("Scan # ", xfile@acquisitionNum[best.peak.scan.num], " @ ",
                   round(xfile@scantime[best.peak.scan.num]/60,1)," min; NL:",
                   formatC(int.max, digits=2,format="e"), sep = ""))
    } else {
      plot(0,0,xlab="",ylab="",main=paste("No quantified peak(s)") )
      plot(0,0,xlab="",ylab="",main=paste("Empty ms1 spectrum") )
      #}
    }
#--------------------------------------------------------medium/heavy--------------------------------------------------------

    NP.value[j] <- paste(n_ms2_peaks, n.candidate.peaks,
                         format(max(k.ms1.int.light.v), digits=1, scientific=T),
                         format(noise.light, digits=1, scientific=T),
						 format(max(k_ms1_int_medium_v), digits=1, scientific=T),
                         format(noise.medium, digits=1, scientific=T),
                         format(max(k_ms1_int_heavy_v), digits=1, scientific=T),
                         format(noise.heavy, digits=1, scientific=T),
                         sep="/")
  } ## each ratio j
  tt <- paste("Entry ", as.character(i), "-  Charge: ", as.character(charge),
              " - M/Z: ", as.character(format(mz.light, digits=7)),
              "and", as.character(format(mz.heavy,digits=7)))
  mtext(tt, line=3.5, outer=T)
  mtext(paste(cross.table[i,"peptide"],"; Mono.mass: ", as.character(mono.mass), "; Mono.mz: ", as.character(round(mono.mz,5)),sep=""),
        cex=0.8, line=2, outer=T)
  mtext(paste(cross.table[i,"ipi"],description),line=0.8, cex=0.8,out=T)
  ## save data in outdf
  lnk.i <- ceiling(i/500)-1
  lnk.j <- (i-1)%%500
  lnk.name <- paste('./PNG/', lnk.i, '/', out.filename.base,'.', lnk.i, '-', lnk.j,'.png',sep='')
  this.df <- c(i, ipi, description, symbol, peptide, round(mass,digits=4), charge, segment,
               i_ratios_LH,i_ratios_LM,i_ratios_HM,l_ratios_LH,l_ratios_LM ,l_ratios_HM,light.int.v, NP.value, r2.v_LH,r2.v_LM,r2.v_HM, entry.index,
               paste('=HYPERLINK(\"./PNG/', lnk.i, '/', out.filename.base,'.', lnk.i, '_', lnk.j,'.png\")',sep=''))
  names(this.df) <- column.names
  out.df <- rbind(out.df, this.df)
  dev.off()
} ## each entry i
#dev.off()

all_table <- out.df
all_table_out <- all_table[F,]
rsq.cutoff <- r2.cutoff

## go from high concentration to low concentration,
## first apply R2 cutoff and sort by IR values
for ( s in seq(ncross, 1) ) {
  colname.R2.LH <- linear.regression.R2.LH[s]
  colname.R2.LM <- linear.regression.R2.LM[s]
  colname.R2.HM <- linear.regression.R2.HM[s]
  colname.IR.LH <- integrated.area.ratio.LH[s]
  colname.IR.LM <- integrated.area.ratio.LM[s]
  colname.IR.HM <- integrated.area.ratio.HM[s]
  rsq.filter <- (all_table[,colname.R2.LH] >= rsq.cutoff & !is.na(all_table[,colname.R2.LH])) | (all_table[,colname.R2.LM] >= rsq.cutoff & !is.na(all_table[,colname.R2.LM])) |(all_table[,colname.R2.HM] >= rsq.cutoff & !is.na(all_table[,colname.R2.HM]))
  table <- all_table[rsq.filter,]
  if (is.vector(table)) {
    table <- data.frame(as.list(table))
  }
  s1 <- as.numeric(table[,colname.IR.LH])
  s2 <- as.numeric(table[,colname.IR.LM])
  s3 <- as.numeric(table[,colname.IR.HM])
  s4 <- as.numeric(table[,"entry"])
  s5 <- as.numeric(table[,"charge"])
  s6 <- as.numeric(table[,"segment"])
  ii <- order(s1, s2, s3, s4, s5, s6)
  table <- table[ii,]
  all_table_out <- rbind(all_table_out, table)
  all_table <- all_table[!rsq.filter,]
}
all_table_out <- rbind(all_table_out, all_table)
## output the final table
row.names(all_table_out) <- as.character(seq(1:dim(all_table_out)[1]) )
all_table_out[,"index"] <- seq(1:dim(all_table_out)[1])
write.table(all_table_out,file=paste(output.path,"/",out.filename.base,".to_excel.txt",sep=""),
            quote=F, sep="\t", row.names=F,na="0.00")
