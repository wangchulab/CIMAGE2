column.names <- c("ipi","description", "symbol", "sequence", "mass", "charge", "segment","run","scan","HL","key","pep_score")
out.df <- matrix(nrow=0,ncol=length(column.names))

args <- commandArgs(trailingOnly=T)
for (filename in args) {
  cat(paste(filename,"\n",sep=""))

  table <- read.table(filename,header=T, sep=",",comment.char="",stringsAsFactors=F)
  for (i in 1:nrow(table)) {
    full.seq <- paste(table[i,"pep_res_before"],table[i,"pep_seq"],table[i,"pep_res_after"],sep=".")
    mod.stat <- table[i,"pep_var_mod_pos"]
    full.seq.vec <- unlist(strsplit(full.seq,""))
    mod.stat.vec <- unlist(strsplit(mod.stat,""))
    mod.tag <- "none"
    heavy.mod.pos <- light.mod.pos <- integer(0)
    if( length(full.seq.vec) == length(mod.stat.vec)) {
      heavy.mod.pos <- which(mod.stat.vec=='2')
      light.mod.pos <- which(mod.stat.vec=='3')
      if (length(heavy.mod.pos) > 0) {
        mod.tag <- 'heavy'
        full.seq.vec[heavy.mod.pos] <- "C*"
      }
      if (length(light.mod.pos) > 0) {
        mod.tag <- 'light'
        full.seq.vec[light.mod.pos] <- "C*"
      }
    }
    full.seq <- paste(full.seq.vec,sep="",collapse="")
    if (mod.tag == "none") next
    scan.title <- unlist(strsplit(as.character(table[i,"pep_scan_title"])," "))[1]
    scan.title.vec <- unlist(strsplit(scan.title,".",fixed=T))
    scan.num <- scan.title.vec[2]
    filename <- scan.title.vec[1]
    filename.vec <- unlist(strsplit(filename,"_"))
    nv <- length(filename.vec)
    segment <- filename.vec[nv]
    run.name <- paste(filename.vec[1:(nv-1)],sep="",collapse="_")
    ipi <- table[i,"prot_acc"]
    full.description <- table[i,"prot_desc"]
    mass <- table[i,"pep_calc_mr"]
    charge <- table[i, "pep_exp_z"]
    pep.score <- table[i,"pep_score"]
    description <- unlist(strsplit(full.description,"OS="))[1]
    symbol.tmp <- unlist(strsplit(full.description,"GN="))[2]
    symbol <- unlist(strsplit(symbol.tmp," "))[1]
    key <- paste(ipi,full.seq,charge,segment,sep=":")
                                        #column.names <- c("ipi","description", "symbol", "sequence", "mass", "charge", "segment","run","scan","HL","key")
    this.df <- c(ipi, description, symbol, full.seq, mass, charge, segment, run.name, scan.num, mod.tag, key, pep.score)
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

ipi.name.table <- as.matrix(levels(as.factor(paste(out.df[,"ipi"],paste(out.df[,"symbol"],out.df[,"description"],sep=" "),sep="\t"))),ncol=1)
colnames(ipi.name.table) <- "name"
write.table(ipi.name.table,file="ipi_name.table", quote=F,sep="\t",row.names=F,na="0.00")
