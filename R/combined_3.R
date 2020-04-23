uniq.tryptic.sequence <- function( raw.sequence ) {
  ## find the smallest tryptic fragment with the labeling sites
  seq.vec <- unlist(strsplit(raw.sequence,""))
  label.pos <- which(seq.vec=="*")
  if ( length(label.pos) == 0 ) {
    ## no mod site, return center sequence
    return( strsplit(raw.sequence,".",fixed=T)[[1]][2] )
  } else {
    first.label <- min( label.pos )
    start <- first.label - 1
    while(seq.vec[start] != ".") {
      start <- start - 1
      if ( (seq.vec[start] == "R" | seq.vec[start] == "K") ) { ##& seq.vec[start+1] != 'P' ) {
        break
      }
    }
    last.label <- max( label.pos )
    end <- last.label+1
    while( seq.vec[end] != "." ) {
      end <- end + 1
      if ( (seq.vec[end-1] == "R" | seq.vec[end-1] == "K") ) { ##& seq.vec[end] != 'P' ) {
        break
      }
    }
    uniq.seq.vec <- seq.vec[(start+1):(end-1)]
    if (length(uniq.seq.vec) == 3 ) { ## only two residues long, return the whole sequence
      xxseq <- strsplit(raw.sequence,".",fixed=T)[[1]][2]
      xxseq.vec <- unlist(strsplit(xxseq,""))
      if (length(label.pos) == 1) {
        return( paste(xxseq.vec[xxseq.vec != "*"], sep="", collapse=""))
      } else {
        return( xxseq )
      }
    }
    if ( length(label.pos) == 1 ) {
      return( paste( uniq.seq.vec[uniq.seq.vec != "*"],sep="",collapse="") )
    } else {
      return( paste( uniq.seq.vec,sep="",collapse="") )
    }
  }
}

## file name from input args
args <- commandArgs(trailingOnly=T)

input.file <- args[1]
dirs <- args[-1]
table <- as.list(dirs)
r2.cutoff <- 0.7

## read in the first table and figure out headers
tmp.table <- read.table(paste(dirs[1],input.file,sep=""),header=T,sep="\t",quote="",comment.char="")
tmp_names <- names(tmp.table)
## to rename columns
v1 <- which(substr(tmp_names,1,6) == "IR_LH.")
nset <- length(v1)
vn1 <- paste("IR_LH.set_",seq(1,nset),sep="")

v2 <- which(substr(tmp_names,1,6) == "IR_LM.")
vn2 <- paste("IR_LM.set_",seq(1,nset),sep="")

v3 <- which(substr(tmp_names,1,6) == "IR_HM.")
vn3 <- paste("IR_HM.set_",seq(1,nset),sep="")

v4 <- which(substr(tmp_names,1,6) == "LR_LH.")
vn4 <- paste("LR_LH.set_",seq(1,nset),sep="")

v5 <- which(substr(tmp_names,1,6) == "LR_LM.")
vn5 <- paste("LR_LM.set_",seq(1,nset),sep="")

v6 <- which(substr(tmp_names,1,6) == "LR_HM.")
vn6 <- paste("LR_HM.set_",seq(1,nset),sep="")

v7 <- which(substr(tmp_names,1,3) == "NP.")
vn7 <- paste("NP.set_",seq(1,nset),sep="")

v8 <- which(substr(tmp_names,1,6) == "R2_LH.")
vn8 <- paste("R2_LH.set_",seq(1,nset),sep="")

v9 <- which(substr(tmp_names,1,6) == "R2_LM.")
vn9 <- paste("R2_LM.set_",seq(1,nset),sep="")

v10 <- which(substr(tmp_names,1,6) == "R2_HM.")
vn10 <- paste("R2_HM.set_",seq(1,nset),sep="")

v11 <- which(substr(tmp_names,1,4) == "INT.")
vn11 <- paste("INT.set_",seq(1,nset),sep="")

nrun <- length(dirs)

all_table <- NULL
for (i in 1:nrun ) {
  table[[i]] <- read.table(paste(dirs[i],input.file,sep=""),header=T,sep="\t",quote="",as.is=T,comment.char="")
  names(table[[i]])[c(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11)] <- c(vn1,vn2,vn3,vn4,vn5,vn6,vn7,vn8,vn9,vn10,vn11)
  table[[i]][,"sequence"] <- as.character( table[[i]][,"sequence"] )
  table[[i]]$run<-i
  table[[i]]$uniq <-i
  table[[i]]$filter <- 0
  for (ii in 1:length(table[[i]]$ipi) ) {
    sequence <- as.character(table[[i]][ii,"sequence"])
    sequence <- uniq.tryptic.sequence( sequence )
    table[[i]]$uniq[ii]<- sequence
  }
  all_table<-rbind(all_table, table[[i]])
}

## set to NA if not passing r2.cutoff
for( i in 1:length(vn1) ) {
  invalid_LH <- (all_table[[vn8[i]]]<r2.cutoff)
  all_table[[vn1[i]]][invalid_LH] <- NA
  invalid_LM <- (all_table[[vn9[i]]]<r2.cutoff)
  all_table[[vn2[i]]][invalid_LM] <- NA
  invalid_HM <- (all_table[[vn10[i]]]<r2.cutoff)
  all_table[[vn3[i]]][invalid_HM] <- NA
}
## only consider entries with at least two valid ratios out of three concentrations
for( i in 1:nrow(all_table) ) {
  all_table[i,"filter"] <- sum(all_table[i,c(vn1[1:nset],vn2[1:nset],vn3[1:nset])]>0, na.rm=T)
}

sp=" "
count <- 0
link.list <- as.list( levels(as.factor(all_table$uniq) ) )
nuniq <- length(link.list)
out_num_matrix <- matrix(0 ,nrow=nuniq,ncol=9*nset)
colnames(out_num_matrix) <- c( paste("mr_LH.set_",seq(1,nset),sep=""),paste("mr_LM.set_",seq(1,nset),sep=""),paste("mr_HM.set_",seq(1,nset),sep=""), paste("lr_LH.set_",seq(1,nset),sep=""), 
                               paste("lr_LM.set_",seq(1,nset),sep=""),paste("lr_HM.set_",seq(1,nset),sep=""),paste("sd_LH.set_",seq(1,nset),sep=""),paste("sd_LM.set_",seq(1,nset),sep=""),
							   paste("sd_HM.set_",seq(1,nset),sep=""))
char.names <- c("index","ipi", "description", "symbol", "sequence", "mass", "run", "charge", "segment", "link")
out.char.matrix <- matrix(" ",nrow=nuniq,ncol=length(char.names))
colnames(out.char.matrix) <- char.names
for (uniq in levels(as.factor(all_table$uniq) ) ) {
  ##ipi <- strsplit(uniq,":")[[1]][1]
  ##seq <- strsplit(uniq,":")[[1]][2]
  match <- all_table[,"uniq"] == uniq  ##(all_table[,"sequence"]==seq) & (all_table$ipi==ipi)
  sub_table <- all_table[match,]
  s1 <- sub_table[,"ipi"]
  s2 <- sub_table[,"sequence"]
  s5 <- -sub_table[,"filter"]
  s3 <- sub_table[,"charge"]
  s4 <- sub_table[,"segment"]
  s6 <- sub_table[,"run"]

  ii <- order(s6,s1,s2,s5,s3,s4)
  count <- count+1
  link.list[[count]] <- which(match)[ii]

  pass <- sub_table$filter>=1
  out.char.matrix[count,"index"] <- as.character(count)
  out.char.matrix[count,"sequence"] <- as.character(uniq)
 # out.char.matrix[count,"run"] <- paste(levels(as.factor(sub_table[,"run"])),sep="",collapse="")
 # out.char.matrix[count,"charge"] <- paste(levels(as.factor(sub_table[,"charge"])),sep="",collapse="")
 # out.char.matrix[count,"segment"] <- paste(levels(as.factor(sub_table[,"segment"])),sep="",collapse="")
  if (sum(pass)>=1) {
    out.char.matrix[count,"run"] <- paste(levels(as.factor(sub_table[pass,"run"])),sep="",collapse="")
    out.char.matrix[count,"charge"] <- paste(levels(as.factor(sub_table[pass,"charge"])),sep="",collapse="")
    out.char.matrix[count,"segment"] <- paste(levels(as.factor(sub_table[pass,"segment"])),sep="",collapse="")
  } else {
    out.char.matrix[count,"run"] <- paste(levels(as.factor(sub_table[,"run"])),sep="",collapse="")
    out.char.matrix[count,"charge"] <- paste(levels(as.factor(sub_table[,"charge"])),sep="",collapse="")
    out.char.matrix[count,"segment"] <- paste(levels(as.factor(sub_table[,"segment"])),sep="",collapse="")
  }
  for ( k in 1:nset ) {
    kk <- k + nset
	kkk <- k + 2*nset
	kkkk <- k + 3*nset
	kkkkk <- k + 4*nset
	kkkkkk <- k + 5*nset
	kkkkkkk <- k + 6*nset
	kkkkkkkk <- k + 7*nset
	kkkkkkkkk <- k + 8*nset
    if (nrun >1) {
      median.per.run <- rep(0,length=nrun)
	  medianlinear.per.run <- rep(0,length=nrun)
      for (dd in 1:nrun) {
        median_LH.per.run[dd] <- round(median(sub_table[pass&(sub_table[,"run"]==dd),vn1[k]],na.rm=T),digits=2)
		median_LM.per.run[dd] <- round(median(sub_table[pass&(sub_table[,"run"]==dd),vn2[k]],na.rm=T),digits=2)
		medianlinear_LH.per.run[dd] <- round(median(sub_table[pass&(sub_table[,"run"]==dd),vn3[k]],na.rm=T),digits=2)
		medianlinear_LM.per.run[dd] <- round(median(sub_table[pass&(sub_table[,"run"]==dd),vn4[k]],na.rm=T),digits=2)
#        if (median.per.run[dd] == 0) {median.per.run[dd] <- NA}
      }
      nrun.valid <- sum( !is.na(median.per.run))
      out_num_matrix[count,k]  <- round(mean(median.per.run,na.rm=T),digits=2)
      out_num_matrix[count,kk] <- round(mean(medianlinear.per.run,na.rm=T),digits=2)
      out_num_matrix[count,kkk] <- round(sd(median.per.run,na.rm=T)+0.01*(nrun.valid-1),digits=2) ## to differentiate multiple 15 ratios from single replicate
	} else  {
      out_num_matrix[count,k]  <- round(median(sub_table[pass,vn1[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kk]  <- round(median(sub_table[pass,vn2[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkk]  <- round(median(sub_table[pass,vn3[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkk]  <- round(median(sub_table[pass,vn4[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkk]  <- round(median(sub_table[pass,vn5[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkk]  <- round(median(sub_table[pass,vn6[k]],na.rm=T),digits=2)
      out_num_matrix[count,kkkkkkk] <- round(sd(sub_table[pass,vn1[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkk]  <- round(sd(sub_table[pass,vn2[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkkk]  <- round(sd(sub_table[pass,vn3[k]],na.rm=T),digits=2)
    }

  }
}

## order by ratio from last set to first set
z.order <- do.call("order", data.frame(out_num_matrix[,seq(nset,1)]))
m.order <- do.call("order", data.frame(out_num_matrix[,2]))
#z.order <- order(data.frame(out_num_matrix[,seq(nset,1)]), decreasing = T)

## draw venn diagrams of averaged ratios
if (nset <=3 ) {
  library(limma)
  png("combined_vennDiagram.png")
  venn.out.matrix <- ! is.na(out_num_matrix[,1:2])
  vc <- vennCounts(venn.out.matrix)
  vennDiagram(vc,main="Number of peptides with valid ratios",counts.col="red")
  uniq_protein_num <- length(levels(as.factor(all_table$sequence) ))
  text(2.4,2.3, labels= uniq_protein_num,cex = 1.3, col="red")
  text(2.4,2, labels= round((vc[2,3] +vc[3,3]+vc[4,3])/sum(vc[,3]),2),cex=1.3, col="red")
  dev.off()
}
##

new_num_matrix <- out_num_matrix[F,]
new.char.matrix <- out.char.matrix[F,]

for ( m in 1:length(z.order)) {
  ii <- z.order[m]
  index <- as.numeric(out.char.matrix[ii,"index"])
  match <- link.list[[index]]
  sub_table <- all_table[match,]

  links <- sub_table[,"link"]
  runs <- as.numeric(sub_table[,"run"])

  this_n_entry <- out_num_matrix[ii,]
  new_num_matrix <- rbind(new_num_matrix,this_n_entry)
  this.c.entry <- out.char.matrix[ii,]
  this.c.entry["index"] <- m
  new.char.matrix <- rbind(new.char.matrix,this.c.entry)

  for ( l in 1:length(links) ) {
    linkfile <- strsplit(as.character(links[l]),'"')[[1]]
    new.filename <- paste(dirs[runs[l]],linkfile[2],sep="")
    new.count <- paste(m, l, sep=".")
    new.link <- paste('=HYPERLINK(\"',new.filename,'\",\"',new.count,'\")',sep='')
    ## fill in information from subtable
    for ( c in char.names ) {
      this.c.entry[c] <- as.character(sub_table[l,c])
    }
    this.c.entry["index"] <- sp
    this.c.entry["link"] <- new.link
    new.char.matrix <- rbind(new.char.matrix,this.c.entry)
    #for ( n in 1:nset ) {
      this_n_entry[1] <- sub_table[l,vn1]
      this_n_entry[2] <- sub_table[l,vn2]
	  this_n_entry[3] <- sub_table[l,vn3]
      this_n_entry[4] <- sub_table[l,vn4]
	  this_n_entry[5] <- sub_table[l,vn5]
	  this_n_entry[6] <- sub_table[l,vn6]
	  this_n_entry[7] <- sub_table[l,vn7]
      this_n_entry[8] <- sub_table[l,vn11]
	  #  ##this_n_entry[n+nset] <- NA
    #}
    new_num_matrix <- rbind(new_num_matrix,this_n_entry)
  }
}

## insert ratio mean and sd in between "mass" and "run" columns
cmass <- which(char.names=="mass")
html.table <- cbind(new.char.matrix[,seq(1,cmass)], ##count to mass
                    new_num_matrix, ## mr and sd
                    new.char.matrix[,seq(cmass+1, length(char.names))] ## run to link
                    )

write.table(html.table,file="combined.txt", quote=F, sep="\t", row.names=F,na="0.00")

png("combined_histogram_IR_LH.png")
ratio <- out_num_matrix[z.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
# linear regression ratio is at the 1st column
 ratio <- ratio[valid,seq(1,nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
# find & remove outliers
#outliers <- boxplot(ratio)$out
#outliers_percent <- paste(length(outliers),length(ratio_original), sep="/")
#ratio <- setdiff(ratio, outliers)
# calculate the mean and sd
mu <- mean(ratio)
sigma <- sd(ratio)
reference.data <- rnorm(length(ratio),mu, sigma)
r.order <- do.call("order", data.frame(reference.data))
reference.data.new <- reference.data[r.order]
x.lm <- lsfit( x=ratio, y=reference.data,intercept=F )
r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
#p_value <- ks.test(reference.data, ratio)$p.value
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
#simulate gaussian distribution and add to the picture
x = seq(from=0,to=2.0,by=0.005)
pdf = dnorm(x, mu, sigma)
lines(pdf ~ x,xlim=c(0,2), col = "green")
#add vertical line at the mean position
abline(v= mu,col="red")
#add comment for the mean,sd and R2
text(1.5,2, labels= paste("mean=", round(mu,3), ";\nsd=", round(sigma,3), ";\nR2=", r2,sep=""))
#title("combined_histogram_IR")
dev.off()

png("combined_histogram_IR_LM.png")
ratio <- out_num_matrix[m.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset])
}
# linear regression ratio is at the 1st column
 ratio <- ratio[valid,seq(1+nset,nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
# find & remove outliers
#outliers <- boxplot(ratio)$out
#outliers_percent <- paste(length(outliers),length(ratio_original), sep="/")
#ratio <- setdiff(ratio, outliers)
# calculate the mean and sd
mu <- mean(ratio)
sigma <- sd(ratio)
reference.data <- rnorm(length(ratio),mu, sigma)
r.order <- do.call("order", data.frame(reference.data))
reference.data.new <- reference.data[r.order]
x.lm <- lsfit( x=ratio, y=reference.data,intercept=F )
r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
#p_value <- ks.test(reference.data, ratio)$p.value
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
#simulate gaussian distribution and add to the picture
x = seq(from=0,to=2.0,by=0.005)
pdf = dnorm(x, mu, sigma)
lines(pdf ~ x,xlim=c(0,2), col = "green")
#add vertical line at the mean position
abline(v= mu,col="red")
#add comment for the mean,sd and R2
text(1.5,2, labels= paste("mean=", round(mu,3), ";\nsd=", round(sigma,3), ";\nR2=", r2,sep=""))
#title("combined_histogram_IR")
dev.off()

png("combined_histogram_LR_LH.png")
ratio <- out_num_matrix[z.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset+nset+nset])
}
# linear regression ratio is at the 1+nest column
ratio <- ratio[valid,seq(1+nset+nset+nset,nset+nset+nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
#calculate the mean and the sd
mu <- mean(ratio)
sigma <- sd(ratio)
#simulate the gaussian distribution data according the mean and sigma
reference.data <- rnorm(length(ratio),mu, sigma)
#re-order the simulated data
r.order <- do.call("order", data.frame(reference.data))
reference.data.new <- reference.data[r.order]
#do linear fitting for the simulated data and observed data
x.lm <- lsfit( x=ratio, y=reference.data.new,intercept=F )
r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
#plot the histogram 
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
#simulate gaussian distribution and add to the picture
x = seq(from=0,to=2.0,by=0.005)
pdf = dnorm(x, mu, sigma)
lines(pdf ~ x,xlim=c(0,2), col = "green")
#add vertical line at the mean position
abline(v= mu,col="red")
#add comment for the mean,sd and R2
text(1.5,2, labels= paste("mean=", round(mu,3), ";\nsd=", round(sigma,3), ";\nR2=", r2, sep=""))
#title("combined_histogram_LR")
dev.off()

png("combined_histogram_LR_LH.png")
ratio <- out_num_matrix[m.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset+nset+nset+nset])
}
# linear regression ratio is at the 1+nest column
ratio <- ratio[valid,seq(1+nset+nset+nset+nset,nset+nset+nset+nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
#calculate the mean and the sd
mu <- mean(ratio)
sigma <- sd(ratio)
#simulate the gaussian distribution data according the mean and sigma
reference.data <- rnorm(length(ratio),mu, sigma)
#re-order the simulated data
r.order <- do.call("order", data.frame(reference.data))
reference.data.new <- reference.data[r.order]
#do linear fitting for the simulated data and observed data
x.lm <- lsfit( x=ratio, y=reference.data.new,intercept=F )
r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
#plot the histogram 
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
#simulate gaussian distribution and add to the picture
x = seq(from=0,to=2.0,by=0.005)
pdf = dnorm(x, mu, sigma)
lines(pdf ~ x,xlim=c(0,2), col = "green")
#add vertical line at the mean position
abline(v= mu,col="red")
#add comment for the mean,sd and R2
text(1.5,2, labels= paste("mean=", round(mu,3), ";\nsd=", round(sigma,3), ";\nR2=", r2, sep=""))
#title("combined_histogram_LR")
dev.off()

png("combined_IR_LH.png")
ratio <- out_num_matrix[z.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1,nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}

x<- seq(nrow(ratio),1)
yl <- c(-4,4) #c(0, max(ratio))
for ( i in 1:nset) {
  plot(x,log2(ratio[,i]),ylim=yl,xlab="Peptide Count",ylab="Observed Ratio(Log2)",col=palette()[i])
  par(new=T)
}
par(new=F)
legend(nrow(ratio)*0.75, yl[2], colnames(ratio), col=palette()[1:nset],pch=1,, text.col=palette()[1:nset])
title("Observed Ratios")
dev.off()

png("combined_IR_LM.png")
ratio <- out_num_matrix[m.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset])
}
ratio <- ratio[valid,seq(1+nset,nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[2]
}

x<- seq(nrow(ratio),1)
yl <- c(-4,4) #c(0, max(ratio))
for ( i in 1:nset) {
  plot(x,log2(ratio[,i]),ylim=yl,xlab="Peptide Count",ylab="Observed Ratio(Log2)",col=palette()[i])
  par(new=T)
}
par(new=F)
legend(nrow(ratio)*0.75, yl[2], colnames(ratio), col=palette()[1:nset],pch=1,, text.col=palette()[1:nset])
title("Observed Ratios")
dev.off()

png("combined_LR_LH.png")
ratio <- out_num_matrix[z.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset+nset+nset])
}
ratio <- ratio[valid,seq(1+nset+nset+nset,nset+nset+nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[3]
}

x<- seq(nrow(ratio),1)
yl <- c(-4,4) #c(0, max(ratio))
for ( i in 1:nset) {
  plot(x,log2(ratio[,i]),ylim=yl,xlab="Peptide Count",ylab="Observed Ratio(Log2)",col=palette()[i])
  par(new=T)
}
par(new=F)
legend(nrow(ratio)*0.75, yl[2], colnames(ratio), col=palette()[1:nset],pch=1,, text.col=palette()[1:nset])
title("Observed Ratios")
dev.off()

png("combined_LR_LM.png")
ratio <- out_num_matrix[m.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset+nset+nset+nset])
}
ratio <- ratio[valid,seq(1+nset+nset+nset+nset,nset+nset+nset+nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[4]
}

x<- seq(nrow(ratio),1)
yl <- c(-4,4) #c(0, max(ratio))
for ( i in 1:nset) {
  plot(x,log2(ratio[,i]),ylim=yl,xlab="Peptide Count",ylab="Observed Ratio(Log2)",col=palette()[i])
  par(new=T)
}
par(new=F)
legend(nrow(ratio)*0.75, yl[2], colnames(ratio), col=palette()[1:nset],pch=1,, text.col=palette()[1:nset])
title("Observed Ratios")
dev.off()
