#segment:in which fraction; run: 150430_3t3_GlcNAl_BA_Tev; scan:ms2 scannum; 
column.names <- c("ipi", "sequence", "mass", "charge", "segment","run","scan","HL","key","pep_score")
out.df <- matrix(nrow=0,ncol=length(column.names))

args <- commandArgs(trailingOnly=T)
cat(args)
mod_symbol <- args[1]
args_len <- length(args)
cat(args[2:args_len])
#mascotfiles <- c()
for (filename in args[2:args_len]) {
  cat(paste(filename,"\n",sep=""))

  table <- read.table(filename,header=T, sep="\t",comment.char="",stringsAsFactors=F)
  if (grepl('light', filename)) {
	mod.tag <- 'light'
  } else {
	mod.tag <- 'heavy'
  }
  for (i in 1:nrow(table)) {
    full.seq <- paste('A',table[i,"Modified.sequence"],'B',sep=".")
	#if (grepl("con",ipi)) next
	full.seq <- gsub(paste(c('\\(',mod_symbol,'\\)'),collapse=''),'*',full.seq)
	full.seq <- gsub('_','',full.seq)
	cat(full.seq)
    #scan.title <- unlist(strsplit(as.character(table[i,"pep_scan_title"])," "))[1]
    #scan.title.vec <- unlist(strsplit(scan.title,".",fixed=T))
    scan.num <- table[i,"MS.MS.scan.number"]
    filename <- table[i,"Raw.file"]
    filename.vec <- unlist(strsplit(filename,"_"))
    nv <- length(filename.vec)
    segment <- filename.vec[nv]
    run.name <- paste(filename.vec[1:(nv-1)],sep="",collapse="_")
    ipi <- table[i,"Leading.proteins"]
    if (grepl("CON",ipi)) next
    #full.description <- table[i,"Protein names"]
    mass <- table[i,"MS.MS.m.z"] * table[i,"Charge"]
    charge <- table[i, "Charge"]
    pep.score <- table[i,"Score"]
    #description <- table[i,"Leading.proteins"]
    #symbol <- table[i,"Leading.proteins"]
    key <- paste(ipi,full.seq,charge,segment,sep=":")
                                        #column.names <- c("ipi","description", "symbol", "sequence", "mass", "charge", "segment","run","scan","HL","key")
    this.df <- c(ipi, full.seq, mass, charge, segment, run.name, scan.num, mod.tag, key, pep.score)
    names(this.df) <- column.names
    out.df <- rbind(out.df,this.df)
  }
}
# output 3 tables

run.name <- levels(as.factor(out.df[,"run"]))
cross.scan.table <- matrix(nrow=0,ncol=3)
cross.scan.table2 <- out.df[,c("key","mass","scan","pep_score")]
for (pep.key in levels(as.factor(cross.scan.table2[,"key"]))){
  key.match <- (cross.scan.table2[,"key"]==pep.key)
  if (sum(key.match) == 1) {
    this.entry <- cross.scan.table2[key.match,c("key","mass","scan")]
  } else {
    entries <- cross.scan.table2[key.match,]
    best.i <- order(entries[,"pep_score"])[1]
    this.entry <- entries[best.i,c("key","mass","scan")]
  }
  cross.scan.table <- rbind(cross.scan.table,this.entry)
}

colnames(cross.scan.table) <- c("key","mass",run.name)
write.table(cross.scan.table,file="cross_scan.table", quote=F,sep="\t",row.names=F,na="0.00")

all.scan.table <- out.df[,c("key","run","scan","HL")]
write.table(all.scan.table,file="all_scan.table", quote=F,sep="\t",row.names=F,na="0.00")

ipi.name.table <- as.matrix(levels(as.factor(out.df[,"ipi"])),ncol=1)
colnames(ipi.name.table) <- "name"
write.table(ipi.name.table,file="ipi_name.table", quote=F,sep="\t",row.names=F,na="0.00")
