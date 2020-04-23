uniq_tryptic_sequence <- function( raw.sequence ) {
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
tmp_table <- read.table(paste(dirs[1],input.file,sep=""),header=T,sep="\t",quote="",comment.char="")
tmp_names <- names(tmp_table)
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
  all_table<-rbind(all_table, table[[i]])
}
all_table$uniq <- all_table$ipi
if (use.mylist) {
  entries <- mylist[,1]
} else {
  entries <- all_table$uniq
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
link_list <- as.list( levels(as.factor(entries) ) )
nuniq <- length(link_list)
out_num_matrix <- matrix(0 ,nrow=nuniq,ncol=13*nset)
colnames(out_num_matrix) <- c( paste("mr_LH.set_",seq(1,nset),sep=""), paste("mr_LM.set_",seq(1,nset),sep=""), paste("mr_HM.set_",seq(1,nset),sep=""),paste("lr_LH.set_",seq(1,nset),sep="") , 
                               paste("lr_LM.set_",seq(1,nset),sep=""), paste("lr_HM.set_",seq(1,nset),sep=""),paste("sd_LH.set_",seq(1,nset),sep=""),paste("sd_LM.set_",seq(1,nset),sep=""),
							   paste("sd_GM.set_",seq(1,nset),sep=""),paste("mean_LH.set_",seq(1,nset),sep=""),paste("mean_LM.set_",seq(1,nset),sep=""), paste("mean_HM.set_",seq(1,nset),sep=""),paste("noqp.set_",seq(1,nset),sep=""))
char_names <- c("index","ipi", "description", "symbol", "sequence", "mass", "run", "charge", "segment", "link")
out_char_matrix <- matrix(" ",nrow=nuniq,ncol=length(char_names))
colnames(out_char_matrix) <- char_names
for (uniq in levels(as.factor(entries) ) ) {
  ##ipi <- strsplit(uniq,":")[[1]][1]
  ##seq <- strsplit(uniq,":")[[1]][2]
  match <- all_table[,"uniq"] == uniq  ##(all_table[,"sequence"]==seq) & (all_table$ipi==ipi)
  if (sum(match) == 0 ) next
  sub_table <- all_table[match,]
  s1 <- sub_table[,"ipi"]
  s2 <- sub_table[,"sequence"]
  s5 <- -sub_table[,"filter"]
  s3 <- sub_table[,"charge"]
  s4 <- sub_table[,"segment"]
  s6 <- sub_table[,"run"]
  
  s2_has_methionine <- rep(0,length(s2))
  if (exclude.methionine) {
    for (ss  in 1:length(s3) ) {
      s2_has_methionine[ss] <- has.methionine(s2[ss])
    }
  }    
  ii <- order(s6,s5,s2,s1,s3,s4)
  count <- count+1
  link_list[[count]] <- which(match)[ii]

  pass <- sub_table$filter>=1
  out_char_matrix[count,"index"] <- as.character(count)
  out_char_matrix[count,"ipi"] <- as.character(uniq)
  out_char_matrix[count,"description"] <- sub_table[1,"description"]
  out_char_matrix[count,"symbol"] <- sub_table[1,"symbol"]
  if (sum(pass)>=1) {
    out_char_matrix[count,"run"] <- paste(levels(as.factor(sub_table[pass,"run"])),sep="",collapse="")
    out_char_matrix[count,"charge"] <- paste(levels(as.factor(sub_table[pass,"charge"])),sep="",collapse="")
    out_char_matrix[count,"segment"] <- paste(levels(as.factor(sub_table[pass,"segment"])),sep="",collapse="")
  } else {
    out_char_matrix[count,"run"] <- paste(levels(as.factor(sub_table[,"run"])),sep="",collapse="")
    out_char_matrix[count,"charge"] <- paste(levels(as.factor(sub_table[,"charge"])),sep="",collapse="")
    out_char_matrix[count,"segment"] <- paste(levels(as.factor(sub_table[,"segment"])),sep="",collapse="")
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
	kkkkkkkkkk <- k + 9*nset
	kkkkkkkkkkk <- k + 10*nset
	kkkkkkkkkkkk <- k + 11*nset
	kkkkkkkkkkkkk <- k + 12*nset
    ## for bmartin, exclude singleton cases for median and sd calculation, also remove c-terminal peptides without mass shift.
    #pass <- pass & ! ( (sub_table[,vn1[k]] == 1) & (sub_table[,vn6[k]] == 1) ) & (s2_has_methionine==0)
    if (exclude.singleton) {
      pass2 <- pass & (sub_table[,vn1[k]] != 15) 
    } else {
      pass2 <- pass  
    }
    if (sum(pass2) >= 1 ) {
      if (nrun>1) {
        median_per_run_LH <- rep(0,length=nrun)
		medianlinear_per_run_LM <- rep(0,length=nrun)
        for (dd in 1:nrun) {
          median_per_run_LH[dd] <- round(median(sub_table[pass2&(sub_table[,"run"]==dd),vn1[k]],na.rm=T),digits=2)
		  median_per_run_LM[dd] <- round(median(sub_table[pass2&(sub_table[,"run"]==dd),vn2[k]],na.rm=T),digits=2)
          medianlinear_per_run_LH[dd] <- round(median(sub_table[pass&(sub_table[,"run"]==dd),vn3[k]],na.rm=T),digits=2)
		  medianlinear_per_run_LM[dd] <- round(median(sub_table[pass&(sub_table[,"run"]==dd),vn4[k]],na.rm=T),digits=2)
		}
        nrun.valid <- sum(!is.na(median_per_run))
        # if there are multiple runs comobined, median/sd column is actually mean/sd of medians, mean column is still mean of all peptide ratios
        out_num_matrix[count,k]  <- round(mean(median_per_run_LH,na.rm=T),digits=2)
		out_num_matrix[count,kk]  <- round(mean(median_per_run_LM,na.rm=T),digits=2)
        out_num_matrix[count,kkk]  <- round(mean(medianlinear_per_run_LH,na.rm=T),digits=2)
		out_num_matrix[count,kkkk]  <- round(mean(medianlinear_per_run_LM,na.rm=T),digits=2)
		out_num_matrix[count,kkkkk] <- round(sd(median_per_run_LH,na.rm=T)+0.01*(nrun.valid-1),digits=2)
		out_num_matrix[count,kkkkkk] <- round(sd(median_per_run_LM,na.rm=T)+0.01*(nrun.valid-1),digits=2)
	  } else {
		#ratio_test[count,1] <- paste(sub_table[pass2,vn1[k]],collapse =";")
		#ratio_test[count,2] <- round(median(sub_table[pass2,vn1[k]],na.rm=T),digits=2)
        out_num_matrix[count,k]  <- round(median(sub_table[pass2,vn1[k]],na.rm=T),digits=2)
		out_num_matrix[count,kk]  <- round(median(sub_table[pass2,vn2[k]],na.rm=T),digits=2)
		out_num_matrix[count,kkk]  <- round(median(sub_table[pass2,vn3[k]],na.rm=T),digits=2)
        out_num_matrix[count,kkkk]  <- round(median(sub_table[pass2,vn4[k]],na.rm=T),digits=2)
		out_num_matrix[count,kkkkk]  <- round(median(sub_table[pass2,vn5[k]],na.rm=T),digits=2)
		out_num_matrix[count,kkkkkk]  <- round(median(sub_table[pass2,vn6[k]],na.rm=T),digits=2)
		out_num_matrix[count,kkkkkkk] <- round(sd(sub_table[pass2,vn1[k]],na.rm=T),digits=2)
		out_num_matrix[count,kkkkkkkk] <- round(sd(sub_table[pass2,vn2[k]],na.rm=T),digits=2)
		out_num_matrix[count,kkkkkkkkk] <- round(sd(sub_table[pass2,vn3[k]],na.rm=T),digits=2)
      }
      out_num_matrix[count,kkkkkkkkkk]  <- round(mean(sub_table[pass2,vn1[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkkkkk]  <- round(mean(sub_table[pass2,vn2[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkkkkkk]  <- round(mean(sub_table[pass2,vn3[k]],na.rm=T),digits=2)
      out_num_matrix[count,kkkkkkkkkkkkk] <- sum(pass2)
    } else {
      out_num_matrix[count,k]  <- round(median(sub_table[pass2,vn1[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kk]  <- round(median(sub_table[pass2,vn2[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkk]  <- round(median(sub_table[pass2,vn3[k]],na.rm=T),digits=2)
      out_num_matrix[count,kkkk]  <- round(median(sub_table[pass2,vn4[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkk]  <- round(median(sub_table[pass2,vn5[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkk]  <- round(median(sub_table[pass2,vn6[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkk] <- round(sd(sub_table[pass2,vn1[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkk] <- round(sd(sub_table[pass2,vn2[k]],na.rm=T),digits=2)
      out_num_matrix[count,kkkkkkkkk] <- round(sd(sub_table[pass2,vn3[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkkkk]  <- round(mean(sub_table[pass2,vn1[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkkkkk]  <- round(mean(sub_table[pass2,vn2[k]],na.rm=T),digits=2)
	  out_num_matrix[count,kkkkkkkkkkkk]  <- round(mean(sub_table[pass2,vn3[k]],na.rm=T),digits=2)
    }
  }
}
if (nuniq > count) {
  out_num_matrix <- out_num_matrix[-seq(count+1,nuniq),]
  out_char_matrix <- out_char_matrix[-seq(count+1,nuniq),]
}
## order by ratio from last set to first set
# z_order <- do.call("order", c(data.frame(out_num_matrix[,seq(nset,1)]))) #del descreasing=descending
#z_order <- order(data.frame(out_num_matrix[,seq(nset,1)]),decreasing=descending)
z_order <- order(data.frame(out_num_matrix[,seq(nset,1)]),decreasing=descending)
m_order <- order(data.frame(out_num_matrix[,2]),decreasing=descending)
h_order <- order(data.frame(out_num_matrix[,3]),decreasing=descending)
## draw venn diagrams of averaged ratios
if (nset <=3 ) {
  library(limma)
  png("combined_vennDiagram.png")
  venn.out.matrix <- ! is.na(out_num_matrix[,1:3])
  vc <- vennCounts(venn.out.matrix)
  vennDiagram(vc,main="Number of proteins with valid ratios",counts.col="red")
  uniq_protein_num <- length(levels(as.factor(all_table$uniq) ))
  text(2.4,2.3, labels= uniq_protein_num,cex = 1.3, col="red")
  text(2.4,2, labels= round((vc[2,3] +vc[3,3]+vc[4,3])/sum(vc[,3]),2),cex=1.3, col="red")
  dev.off()
}
##

new_num_matrix <- out_num_matrix[F,]
new_char_matrix <- out_char_matrix[F,]

for ( m in 1:length(z_order)) {
  ii <- z_order[m]
  index <- as.numeric(out_char_matrix[ii,"index"])
  match <- link_list[[index]]
  sub_table <- all_table[match,]

  links <- sub_table[,"link"]
  runs <- as.numeric(sub_table[,"run"])

  this_n_entry <- out_num_matrix[ii,]
  new_num_matrix <- rbind(new_num_matrix,this_n_entry)
  this_c_entry <- out_char_matrix[ii,]
  this_c_entry["index"] <- m
  new_char_matrix <- rbind(new_char_matrix,this_c_entry)

  for ( l in 1:length(links) ) {
    linkfile <- strsplit(as.character(links[l]),'"')[[1]]
    new.filename <- paste(dirs[runs[l]],linkfile[2],sep="")
    new.count <- paste(m, l, sep=".")
    new.link <- paste('=HYPERLINK(\"',new.filename,'\",\"',new.count,'\")',sep='')
    ## fill in information from subtable
    for ( c in char_names ) {
      this_c_entry[c] <- as.character(sub_table[l,c])
    }
    ## empty these columns as they are same as those in the averaged line
    this_c_entry["index"] <- sp
    this_c_entry["ipi"] <- sp
    this_c_entry["description"] <- sp
    this_c_entry["symbol"] <- sp

    this_c_entry["link"] <- new.link
    new_char_matrix <- rbind(new_char_matrix,this_c_entry)
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
cmass <- which(char_names=="mass")
html.table <- cbind(new_char_matrix[,seq(1,cmass)], ##count to mass
                    new_num_matrix, ## mr and sd
                    new_char_matrix[,seq(cmass+1, length(char_names))] ## run to link
                    )

write.table(html.table,file="combined.txt", quote=F, sep="\t", row.names=F,na="0.00")
html.table2 <- html.table[html.table[,"index"]!=sp,]
write.table(html.table2,file="combined_averaged_ratios.txt", quote=F, sep="\t", row.names=F,na="0.00")

png("combined_histogram_IR_LH.png")
ratio <- out_num_matrix[z_order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1,nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_histogram_IR_LM.png")
ratio <- out_num_matrix[m_order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset])
}
ratio <- ratio[valid,seq(1+nset,nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_histogram_IR_HM.png")
ratio <- out_num_matrix[z_order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+2*nset])
}
ratio <- ratio[valid,seq(1+2*nset,3*nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_histogram_LR_LH.png")
ratio <- out_num_matrix[z_order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+3*nset])
}
ratio <- ratio[valid,seq(1+3*nset,4*nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_histogram_LR_LM.png")
ratio <- out_num_matrix[m_order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+4*nset])
}
ratio <- ratio[valid,seq(1+4*nset,5*nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_histogram_LR_HM.png")
ratio <- out_num_matrix[m_order,]
valid <- rep(T,nuniq)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+5*nset])
}
ratio <- ratio[valid,seq(1+5*nset,6*nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
}
hist(ratio,xlim=c(0,2),breaks=seq(min(ratio),max(ratio)+0.02,by=0.02),freq=F)
lines(density(ratio),xlim=c(0,2))
dev.off()

png("combined_IR_LH.png")
ratio <- out_num_matrix[z_order,]
valid <- rep(T,count)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1,nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
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

png("combined_IR_LM.png")
ratio <- out_num_matrix[m_order,]
valid <- rep(T,count)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+nset])
}
ratio <- ratio[valid,seq(1+nset,2*nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
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

png("combined_IR_HM.png")
ratio <- out_num_matrix[h_order,]
valid <- rep(T,count)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i+2*nset])
}
ratio <- ratio[valid,seq(1+2*nset,3*nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
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


png("combined_LR_LH.png")
ratio <- out_num_matrix[z_order,]
valid <- rep(T,count)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1+nset,nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
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

png("combined_LR_LM.png")
ratio <- out_num_matrix[m_order,]
valid <- rep(T,count)
for ( i in 1:nset) {
  valid <- valid & !is.na(ratio[,i])
}
ratio <- ratio[valid,seq(1+nset,nset+nset)]

if ( is.vector(ratio) ) {
  ratio <- matrix( ratio, byrow=T,ncol=1 )
  colnames(ratio) <- colnames(out_num_matrix)[1]
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
