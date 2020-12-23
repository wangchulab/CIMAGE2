## file name from input args
args <- commandArgs(trailingOnly=T)

nset <- length(args)/3
ratio.files <- args[seq(1,by=3,length=nset)]
input.cols  <- args[seq(2,by=3,length=nset)]
output.cols <- args[seq(3,by=3,length=nset)]

table <- as.list(ratio.files)

all.uniq <- NULL
sp=" "
for (i in 1:nset ) {
  tmp.table <- read.table(ratio.files[i], header=T,sep="\t",quote="",as.is=T,comment.char="")
  table[[i]] <- tmp.table[! is.na(tmp.table[,"index"]),]

  table[[i]]$uniq <-i
  for (ii in 1:length(table[[i]]$ipi) ) {
    ipi <- as.character(table[[i]][ii,"ipi"])
    description <- as.character(table[[i]][ii,"description"])
    symbol <- as.character(table[[i]][ii,"symbol"])
    table[[i]]$uniq[ii]<- paste(ipi,description,symbol,sep=":")
  }
  all.uniq<-c(all.uniq,table[[i]]$uniq)
}

count <- 0
link.list <- as.list( levels(as.factor(all.uniq) ) )
nuniq <- length(link.list)
out.num.matrix <- matrix(NA, nrow=nuniq,ncol=4*nset)
colnames(out.num.matrix) <- c(paste(output.cols,"median",sep=""),paste(output.cols,"mean",sep="."),paste(output.cols,"sd",sep="."),paste(output.cols,"noqp",sep="."))
char.names <- c("index","ipi", "description", "symbol", "sequence")
out.char.matrix <- matrix(" ",nrow=nuniq,ncol=length(char.names))
colnames(out.char.matrix) <- char.names
for (uniq in levels(as.factor(all.uniq) ) ) {
  count <- count + 1
  tmp.split <- strsplit(uniq,":")[[1]]
  out.char.matrix[count,"index"] <- as.character(count)
  out.char.matrix[count,"ipi"] <- tmp.split[1]
  out.char.matrix[count,"description"] <- tmp.split[2]
  out.char.matrix[count,"symbol"] <- tmp.split[3]

  for ( i in 1:nset ) {
    match <- table[[i]][,"uniq"] == uniq
    if ( sum(match) == 1 ) {
      ratio <- table[[i]][match,paste("mr",input.cols[i],sep=".")]
      ratio.mean <- table[[i]][match,paste("mean",input.cols[i],sep=".")]
      noqp <- table[[i]][match,paste("noqp",input.cols[i],sep=".")]
      sd <- table[[i]][match,paste("sd",input.cols[i],sep=".")]
      if (ratio > 0.0) {
        out.num.matrix[count,i] <- ratio
        out.num.matrix[count,i+nset] <- ratio.mean
        out.num.matrix[count,i+2*nset] <- sd
        out.num.matrix[count,i+3*nset] <- noqp
      } else {
        out.num.matrix[count,i] <- NA
        out.num.matrix[count,i+nset] <- NA
        out.num.matrix[count,i+2*nset] <- NA
        out.num.matrix[count,i+3*nset] <- NA
      }
    }
  }
}

## order by ratio from last set to first set
z.order <- do.call("order",c(data.frame(out.num.matrix[,seq(1,2*nset)]), na.last=T))

html.table <- cbind(out.char.matrix,out.num.matrix)[z.order,]

html.table[,"index"] <- seq(1,nrow(html.table))

write.table(html.table, file=paste("compare_averaged_ratios",paste(output.cols,sep="",collapse="_"),
                          "txt",sep="."),
            quote=F, sep="\t", row.names=F,na="0.00")

library(limma)
png(paste("compare_averaged_ratios",paste(output.cols,sep="",collapse="_"),"vennDiagram.png",sep="."))
venn.out.matrix <- ! is.na(out.num.matrix[,1:nset])
vc <- vennCounts(venn.out.matrix)
vennDiagram(vc,main="Number of peptides with valid ratios",counts.col="red")
dev.off()
