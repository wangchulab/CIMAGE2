## file name from input args
args <- commandArgs(trailingOnly=T)

nset <- length(args)/3
cat(nset)
ratio.files <- args[seq(1,by=3,length=nset)]
input.cols  <- args[seq(2,by=3,length=nset)]
output.cols <- args[seq(3,by=3,length=nset)]

tmp.table <- table <- as.list(ratio.files)

all.uniq <- NULL
sp=" "
for (i in 1:nset ) {
  tmp.table[[i]] <- read.table(ratio.files[i], header=T,sep="\t",quote="",as.is=T,comment.char="")
  for (j in 1:nrow(tmp.table[[i]])) {
    if ( (! is.na(tmp.table[[i]][j,"index"])) & (tmp.table[[i]][j,"ipi"]==" ") ) {
      tmp.table[[i]][j,"ipi"] = tmp.table[[i]][j+1,"ipi"]
      tmp.table[[i]][j,"description"] = tmp.table[[i]][j+1,"description"]
      tmp.table[[i]][j,"symbol"] = tmp.table[[i]][j+1,"symbol"]
    }
  }
  table[[i]] <- tmp.table[[i]][! is.na(tmp.table[[i]][,"index"]),]

  table[[i]][,"sequence"] <- as.character( table[[i]][,"sequence"] )
  table[[i]]$uniq <-i
  for (ii in 1:length(table[[i]]$ipi) ) {
    ipi <- as.character(table[[i]][ii,"ipi"])
    description <- as.character(table[[i]][ii,"description"])
    symbol <- as.character(table[[i]][ii,"symbol"])
    sequence <- as.character(table[[i]][ii,"sequence"])
    table[[i]]$uniq[ii]<- paste(ipi,description,symbol,sequence,sep="::")
  }
  all.uniq<-c(all.uniq,table[[i]]$uniq)
}

count <- 0
link.list <- as.list( levels(as.factor(all.uniq) ) )
nuniq <- length(link.list)
out.num.matrix <- matrix(NA, nrow=nuniq,ncol=2*nset)
colnames(out.num.matrix) <- c(paste("mr",output.cols,sep="."),paste("sd",output.cols,sep="."))
num.names <- colnames(out.num.matrix)
char.names <- c("index","ipi", "description", "symbol", "sequence")
out.char.matrix <- matrix(" ",nrow=nuniq,ncol=length(char.names))
colnames(out.char.matrix) <- char.names
for (uniq in levels(as.factor(all.uniq) ) ) {
  count <- count + 1
  tmp.split <- strsplit(uniq,"::")[[1]]
  out.char.matrix[count,"index"] <- as.character(count)
  out.char.matrix[count,"ipi"] <- tmp.split[1]
  out.char.matrix[count,"description"] <- tmp.split[2]
  out.char.matrix[count,"symbol"] <- tmp.split[3]
  out.char.matrix[count,"sequence"] <- tmp.split[4]

  for ( i in 1:nset ) {
    match <- table[[i]][,"uniq"] == uniq
    if ( sum(match) == 1 ) {
      ratio <- table[[i]][match,paste("mr",input.cols[i],sep=".")]
      sd <- table[[i]][match,paste("sd",input.cols[i],sep=".")]
      if (!is.na(ratio) & ratio > 0.0) {
        out.num.matrix[count,i] <- ratio
        out.num.matrix[count,i+nset] <- sd
      } else {
        out.num.matrix[count,i] <- NA
        out.num.matrix[count,i+nset] <- NA
      }
    }
  }
}

## order by ratio from last set to first set
z.order <- do.call("order",c(data.frame(out.num.matrix[,seq(1,2*nset)]), na.last=T))

html.table <- cbind(out.char.matrix,out.num.matrix)[z.order,]

html.table[,"index"] <- seq(1,nrow(html.table))

write.table(html.table, file=paste("compare_averaged_ratios",paste(output.cols,sep="",collapse="_"),
                           "to_excel", "txt",sep="."),
            quote=F, sep="\t", row.names=F,na="0.00")

library(limma)
png(paste("compare_averaged_ratios",paste(output.cols,sep="",collapse="_"),"vennDiagram.png",sep="."))
venn.out.matrix <- ! is.na(out.num.matrix[,1:nset])
if (nset <= 5) { 
	vc <- vennCounts(venn.out.matrix)
	vennDiagram(vc,main="Number of peptides with valid ratios",counts.col="red")
}	
dev.off()

if (ncol(venn.out.matrix)==2) {
  png(paste("compare_averaged_ratios",paste(output.cols,sep="",collapse="_"),"png",sep="."))
  out.num.matrix[is.na(out.num.matrix)] <- 0.0
  qtt <- quantile(as.numeric(out.num.matrix),probs=seq(0,1,0.01))
  limit <- c(qtt[2],qtt[length(qtt)-1])
  plot(out.num.matrix[,1], out.num.matrix[,2],main="averaged ratios comparison",xlab=output.cols[1],ylab=output.cols[2],
       xlim=limit,ylim=limit)
  abline(0,1)

  dev.off()
}


### generate a html ouput with raw spectrum link
link.char.matrix <- matrix(" ",nrow=nuniq,ncol=1)
colnames(link.char.matrix) <- "link"
text.table <- cbind(html.table,link.char.matrix)
out.df <- matrix(nrow=0,ncol=ncol(text.table))
ratio.files.prefix <- ratio.files
for (i in 1:nset) {
  xx <- strsplit(ratio.files[i],"/",fixed=T)
  ratio.files.prefix[i] <- paste(head(xx[[1]],-1),collapse='/',sep='')
}
for( ir in 1: nrow(text.table) ) {
  this.df <- text.table[ir,]
  out.df <- rbind(out.df,this.df)
  this.seq <- text.table[ir,"sequence"]
  this.ipi <- text.table[ir,"ipi"]
  for (i in 1:nset) {
    istart <- which( (tmp.table[[i]][,"ipi"]==this.ipi) & ( tmp.table[[i]][,"sequence"]==this.seq) )
    if (length(istart) ==1 ) {
      inext <- istart + 1
      while( is.na(tmp.table[[i]][inext,"index"]) & (tmp.table[[i]][inext,"ipi"]==this.ipi) ) { # this entry matches
        this.ratio <- tmp.table[[i]][inext,paste("mr",input.cols[i],sep=".")]
        this.sd    <- tmp.table[[i]][inext,paste("sd",input.cols[i],sep=".")]
        if (this.ratio > 0) { # specfic ratio > 0
          this.df[char.names] <- " " #space as default
          this.df["sequence"] <- tmp.table[[i]][inext,"sequence"]
          this.df[num.names] <- NA # NA as default
          this.df[num.names[i]] <- this.ratio
          this.df[num.names[i+nset]] <-this.sd
          this.link <- tmp.table[[i]][inext,"link"]
          xx <- strsplit(this.link,'"',fixed=T)[[1]]
          new.filename <- paste(ratio.files.prefix[i],xx[2],sep="/")
          new.count <- paste(i,xx[4],sep=".")
          new.link <- new.link <- paste('=HYPERLINK(\"',new.filename,'\",\"',new.count,'\")',sep='')
          this.df["link"] <- new.link
          out.df <- rbind(out.df, this.df)
        }
        inext <- inext + 1
        if (inext > nrow(tmp.table[[i]])) break
      } # this entry matches
    } # if
  } #nset
} # ir

outname<-paste("compare_averaged_ratios", paste(output.cols, sep="", collapse="_"), "txt", sep=".")

write.table(out.df, file=outname, quote=F, sep="\t", row.names=F,na="0.00")
write(outname,file="outname")
