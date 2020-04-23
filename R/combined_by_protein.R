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
      return( strsplit(raw.sequence,".",fixed=T)[[1]][2] )
    }
    if ( length(label.pos) == 1 ) {
      return( paste( uniq.seq.vec[uniq.seq.vec != "*"],sep="",collapse="") )
    } else {
      return( paste( uniq.seq.vec,sep="",collapse="") )
    }
  }
}

has.methionine <- function( raw.sequence ) {
  ## find if the sequence contains a methionine or not
  seq.vec <- unlist(strsplit(strsplit(raw.sequence,".",fixed=T)[[1]][2],""))
  label.pos <- which(seq.vec=="M")
  return(length(label.pos)) 
}

## file name from input args
args <- commandArgs(trailingOnly=T)
exclude.singleton <- F
descending <- F
exclude.methionine <- F
mylist <- ""
use.mylist <- F
if ( args[1] == "exclude_singleton" ) {
  exclude.singleton <- T
  args <- args[-1]
}
if ( args[1] == "exclude_MET" ) {
  exclude.methionine <- T
  args <- args[-1]
}
if ( args[1] == "descending" ) {
  descending <- T
  args <- args[-1]
}
if ( args[1] == "mylist_none" ) {
  use.mylist <- F
} else {
  use.mylist <- T
  mylist <- read.table(args[1],header=F)
}
args <- args[-1]

input.file <- args[1]
dirs <- args[-1]
table <- as.list(dirs)
r2.cutoff <- 0.8

## read in the first table and figure out headers
tmp.table <- read.table(paste(dirs[1],input.file,sep=""),header=T,sep="\t",quote="",comment.char="")
tmp.names <- names(tmp.table)
## to rename columns
v1 <- which(substr(tmp.names,1,3) == "IR.")
nset <- length(v1)
vn1 <- paste("IR.set_",seq(1,nset),sep="")
v2 <- which(substr(tmp.names,1,3) == "LR.")
vn2 <- paste("LR.set_",seq(1,nset),sep="")

v3 <- which(substr(tmp.names,1,3) == "NP.")
vn3 <- paste("NP.set_",seq(1,nset),sep="")
v4 <- which(substr(tmp.names,1,3) == "R2.")
vn4 <- paste("R2.set_",seq(1,nset),sep="")
v5 <- which(substr(tmp.names,1,4) == "INT.")
vn5 <- paste("INT.set_",seq(1,nset),sep="")

nrun <- length(dirs)
all.table <- NULL
for (i in 1:nrun ) {
  table[[i]] <- read.table(paste(dirs[i],input.file,sep=""),header=T,sep="\t",quote="",as.is=T,comment.char="")
  names(table[[i]])[c(v1,v2,v3,v4,v5)] <- c(vn1,vn2,vn3,vn4,vn5)
  table[[i]][,"sequence"] <- as.character( table[[i]][,"sequence"] )
  table[[i]]$run<-i
  table[[i]]$uniq <-i
  table[[i]]$filter <- 0
  all.table<-rbind(all.table, table[[i]])
}
all.table$uniq <- all.table$ipi
if (use.mylist) {
  entries <- mylist[,1]
} else {
  entries <- all.table$uniq
}
## set to NA if not passing r2.cutoff
for( i in 1:length(vn1) ) {
  invalid <- (all.table[[vn4[i]]]<r2.cutoff)
  all.table[[vn1[i]]][invalid] <- NA
}
## only consider entries with at least two valid ratios out of three concentrations
for( i in 1:nrow(all.table) ) {
  all.table[i,"filter"] <- sum(all.table[i,vn1[1:nset]]>0, na.rm=T)
}

sp=" "
count <- 0
link.list <- as.list( levels(as.factor(entries) ) )
nuniq <- length(link.list)
out.num.matrix <- matrix(0 ,nrow=nuniq,ncol=5*nset)
colnames(out.num.matrix) <- c( paste("mr.set_",seq(1,nset),sep=""), paste("lr.set_",seq(1,nset),sep="") ,paste("sd.set_",seq(1,nset),sep=""),paste("mean.set_",seq(1,nset),sep=""), paste("noqp.set_",seq(1,nset),sep=""))
char.names <- c("index","ipi", "description", "symbol", "sequence", "mass", "run", "charge", "segment", "link")
out.char.matrix <- matrix(" ",nrow=nuniq,ncol=length(char.names))
colnames(out.char.matrix) <- char.names
for (uniq in levels(as.factor(entries) ) ) {
  ##ipi <- strsplit(uniq,":")[[1]][1]
  ##seq <- strsplit(uniq,":")[[1]][2]
  match <- all.table[,"uniq"] == uniq  ##(all.table[,"sequence"]==seq) & (all.table$ipi==ipi)
  if (sum(match) == 0 ) next
  sub.table <- all.table[match,]
  s1 <- sub.table[,"ipi"]
  s2 <- sub.table[,"sequence"]
  s5 <- -sub.table[,"filter"]
  s3 <- sub.table[,"charge"]
  s4 <- sub.table[,"segment"]
  s6 <- sub.table[,"run"]
  
  s2.has.methionine <- rep(0,length(s2))
  if (exclude.methionine) {
    for (ss  in 1:length(s3) ) {
      s2.has.methionine[ss] <- has.methionine(s2[ss])
    }
  }    
  ii <- order(s6,s5,s2,s1,s3,s4)
  count <- count+1
  link.list[[count]] <- which(match)[ii]

  pass <- sub.table$filter>=1
  out.char.matrix[count,"index"] <- as.character(count)
  out.char.matrix[count,"ipi"] <- as.character(uniq)
  out.char.matrix[count,"description"] <- sub.table[1,"description"]
  out.char.matrix[count,"symbol"] <- sub.table[1,"symbol"]
  if (sum(pass)>=1) {
    out.char.matrix[count,"run"] <- paste(levels(as.factor(sub.table[pass,"run"])),sep="",collapse="")
    out.char.matrix[count,"charge"] <- paste(levels(as.factor(sub.table[pass,"charge"])),sep="",collapse="")
    out.char.matrix[count,"segment"] <- paste(levels(as.factor(sub.table[pass,"segment"])),sep="",collapse="")
  } else {
    out.char.matrix[count,"run"] <- paste(levels(as.factor(sub.table[,"run"])),sep="",collapse="")
    out.char.matrix[count,"charge"] <- paste(levels(as.factor(sub.table[,"charge"])),sep="",collapse="")
    out.char.matrix[count,"segment"] <- paste(levels(as.factor(sub.table[,"segment"])),sep="",collapse="")
  }
  for ( k in 1:nset ) {
    kk <- k + nset
    kkk <- kk + nset
    kkkk <- kkk + nset
	kkkkk <- kkkk + nset
    ## for bmartin, exclude singleton cases for median and sd calculation, also remove c-terminal peptides without mass shift.
    pass <- pass & ! ( (sub.table[,vn1[k]] == 1) & (sub.table[,vn3[k]] == 1) ) & (s2.has.methionine==0)
    if (exclude.singleton) {
      pass2 <- pass & (sub.table[,vn1[k]] != 15) 
    } else {
      pass2 <- pass  
    }
    if (sum(pass2) >= 1 ) {
      if (nrun>1) {
        median.per.run <- rep(0,length=nrun)
		medianlinear.per.run <- rep(0,length=nrun)
        for (dd in 1:nrun) {
          median.per.run[dd] <- round(median(sub.table[pass2&(sub.table[,"run"]==dd),vn1[k]],na.rm=T),digits=2)
          medianlinear.per.run[dd] <- round(median(sub.table[pass&(sub.table[,"run"]==dd),vn2[k]],na.rm=T),digits=2)
		}
        nrun.valid <- sum(!is.na(median.per.run))
        # if there are multiple runs comobined, median/sd column is actually mean/sd of medians, mean column is still mean of all peptide ratios
        out.num.matrix[count,k]  <- round(mean(median.per.run,na.rm=T),digits=2)
        out.num.matrix[count,kk]  <- round(mean(medianlinear.per.run,na.rm=T),digits=2)
		out.num.matrix[count,kkk] <- round(sd(median.per.run,na.rm=T)+0.01*(nrun.valid-1),digits=2)
	  } else {
        out.num.matrix[count,k]  <- round(median(sub.table[pass2,vn1[k]],na.rm=T),digits=2)
        out.num.matrix[count,kk]  <- round(median(sub.table[pass2,vn2[k]],na.rm=T),digits=2)
		out.num.matrix[count,kkk] <- round(sd(sub.table[pass2,vn1[k]],na.rm=T),digits=2)
      }
      out.num.matrix[count,kkkk]  <- round(mean(sub.table[pass2,vn1[k]],na.rm=T),digits=2)
      out.num.matrix[count,kkkkk] <- sum(pass2)
    } else {
      out.num.matrix[count,k]  <- round(median(sub.table[pass,vn1[k]],na.rm=T),digits=2)
      out.num.matrix[count,kk]  <- round(median(sub.table[pass,vn2[k]],na.rm=T),digits=2)
	  out.num.matrix[count,kkk] <- round(sd(sub.table[pass,vn1[k]],na.rm=T),digits=2)
      out.num.matrix[count,kkkk]  <- round(mean(sub.table[pass,vn1[k]],na.rm=T),digits=2)
      out.num.matrix[count,kkkkk] <- sum(pass2)
    }
  }
}
if (nuniq > count) {
  out.num.matrix <- out.num.matrix[-seq(count+1,nuniq),]
  out.char.matrix <- out.char.matrix[-seq(count+1,nuniq),]
}
## order by ratio from last set to first set
# z.order <- do.call("order", c(data.frame(out.num.matrix[,seq(nset,1)]))) #del descreasing=descending
#z.order <- order(data.frame(out.num.matrix[,seq(nset,1)]),decreasing=descending)
z.order <- order(data.frame(out.num.matrix[,seq(nset,1)]),decreasing=descending)
## draw venn diagrams of averaged ratios
if (nset <=3 ) {
  library(limma)
  png("combined_vennDiagram.png")
  venn.out.matrix <- ! is.na(out.num.matrix[,seq(1,length=nset)])
  vc <- vennCounts(venn.out.matrix)
  vennDiagram(vc,main="Number of proteins with valid ratios",counts.col="red")
  uniq_sequence_num <- length(levels(as.factor(all.table$sequence) ))
  text(2.3,2, labels= uniq_sequence_num,cex = 1.3, col="red")
  text(2.3,1.7, labels= round(vc[2,2]/(vc[2,2]+vc[1,2]),2),cex=1.3, col="red")
  dev.off()
}
##

new.num.matrix <- out.num.matrix[F,]
new.char.matrix <- out.char.matrix[F,]

for ( m in 1:length(z.order)) {
  ii <- z.order[m]
  index <- as.numeric(out.char.matrix[ii,"index"])
  match <- link.list[[index]]
  sub.table <- all.table[match,]

  links <- sub.table[,"link"]
  runs <- as.numeric(sub.table[,"run"])

  this.n.entry <- out.num.matrix[ii,]
  new.num.matrix <- rbind(new.num.matrix,this.n.entry)
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
      this.c.entry[c] <- as.character(sub.table[l,c])
    }
    ## empty these columns as they are same as those in the averaged line
    this.c.entry["index"] <- sp
    this.c.entry["ipi"] <- sp
    this.c.entry["description"] <- sp
    this.c.entry["symbol"] <- sp

    this.c.entry["link"] <- new.link
    new.char.matrix <- rbind(new.char.matrix,this.c.entry)
    for ( n in 1:nset ) {
      this.n.entry[n] <- sub.table[l,vn1[n]]
      this.n.entry[n+nset] <- as.character(sub.table[l,vn2[n]])
      ##this.n.entry[n+nset] <- NA
    }
    new.num.matrix <- rbind(new.num.matrix,this.n.entry)
  }
}

## insert ratio mean and sd in between "mass" and "run" columns
cmass <- which(char.names=="mass")
html.table <- cbind(new.char.matrix[,seq(1,cmass)], ##count to mass
                    new.num.matrix, ## mr and sd
                    new.char.matrix[,seq(cmass+1, length(char.names))] ## run to link
                    )

write.table(html.table,file="combined.txt", quote=F, sep="\t", row.names=F,na="0.00")
html.table2 <- html.table[html.table[,"index"]!=sp,]
write.table(html.table2,file="combined_averaged_ratios.txt", quote=F, sep="\t", row.names=F,na="0.00")

png("combined_histogram_IR.png")
ratio <- out.num.matrix[z.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1,nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out.num.matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_histogram_LR.png")
ratio <- out.num.matrix[z.order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1+nset,nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out.num.matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_IR.png")
ratio <- out.num.matrix[z.order,]
valid <- rep(T,count)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1,nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out.num.matrix)[1]
}
x<- seq(nrow(ratio),1)
yl <- c(-4,4)# (c(0, max(ratio))
for ( i in 1:nset) {
  plot(x,log2(ratio[,i]),ylim=yl,xlab="Peptide Count",ylab="Observed Ratio(Log2)",col=palette()[i])
  par(new=T)
}
par(new=F)
legend(nrow(ratio)*0.75, yl[2], colnames(ratio), col=palette()[1:nset],pch=1,, text.col=palette()[1:nset])
title("Observed Ratios")
dev.off()

png("combined_LR.png")
ratio <- out.num.matrix[z.order,]
valid <- rep(T,count)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1+nset,nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out.num.matrix)[1]
}
x<- seq(nrow(ratio),1)
yl <- c(-4,4)# (c(0, max(ratio))
for ( i in 1:nset) {
  plot(x,log2(ratio[,i]),ylim=yl,xlab="Peptide Count",ylab="Observed Ratio(Log2)",col=palette()[i])
  par(new=T)
}
par(new=F)
legend(nrow(ratio)*0.75, yl[2], colnames(ratio), col=palette()[1:nset],pch=1,, text.col=palette()[1:nset])
title("Observed Ratios")
dev.off()
