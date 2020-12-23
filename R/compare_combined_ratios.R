args <- commandArgs(trailingOnly=T)
if (length(args) < 6 ) {
  print("R --args table1 column1 name1 table2 column2 name2")
  quit()
}
table1 <- args[1]
col1 <- args[2]
name1 <- args[3]
table2 <- args[4]
col2 <- args[5]
name2 <- args[6]

cwd <- getwd()
if (substr(table1,1,1) != "/" ) {
  table1 <- paste(cwd,table1,sep="/")
}
if (substr(table1,1,1) != "/" ) {
  table2 <- paste(cwd,table2,sep="/")
}

par(oma=c(0,0,5,0))
##postscript("compare_conc_time.ps")
outname <- paste("compare_ratios",name1,"vs",name2,sep=".")
png(paste(outname,"png",sep="."))

table.conc <- read.table(table1,header=T,sep="\t",stringsAsFactors=F)
table.time <- read.table(table2,header=T,sep="\t",stringsAsFactors=F)

table.conc <- table.conc[(table.conc[,"index"]!=" "), ]
table.time <- table.time[(table.time[,"index"]!=" "), ]

seq.conc <- table.conc[,"sequence"]
seq.time <- table.time[,"sequence"]

seq.all <- c(seq.conc,seq.time)
seq <- character(0)
r.conc <- r.time <- numeric(0)

for ( s in levels(as.factor(seq.all)) ) {
  seq <- c(seq,s)
  i.conc <- which(seq.conc==s)
  if (length(i.conc) == 1 ) {
    r.conc <- c(r.conc,table.conc[i.conc,col1])
  } else {
    r.conc <- c(r.conc,0)
  }
  i.time <- which(seq.time==s)
  if (length(i.time) == 1 ) {
    r.time <- c(r.time,table.time[i.time,col2])
  } else {
    r.time <- c(r.time,0)
  }
}
limit <- range(c(r.conc,r.time))

common <- r.conc>0 & r.time>0
fit.conc <- r.conc[common]
fit.time <- r.time[common]
x.lm <- lsfit( x=fit.conc, y=fit.time, intercept=F)
r2 <- round(as.numeric(ls.print(x.lm,print.it=F)$summary[1,2]),digits=2)
cor <- x.lm$coefficients
plot(r.conc,r.time,xlab=paste("ratio from",name1), ylab=paste("ratio from", name2), xlim=limit, ylim=limit, main=paste(name1,"vs", name2),sub=paste("Xcor=",round(cor,2),"; R2=",round(r2,2)))
abline(0,cor)
dev.off()

out.table <- cbind(fit.conc, fit.time)
colnames(out.table) <- c(name1, name2)
rownames(out.table) <- seq[common]
write.table(out.table, file=paste(outname,"txt",sep="."), sep="\t",row.names=T)
