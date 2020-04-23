read.input.params <- function( file.name ) {
  params <- list()
  raw.input <- scan(file.name, what=character(), quiet=TRUE, comment.char="!")
  equal.pos <- which( raw.input == "=" )
  n.params <- length(equal.pos)
  ## param 1:(N-1)
  for ( i in 1:(n.params-1) ) {
    param.name <- raw.input[equal.pos[i]-1]
    param.data <- raw.input[(equal.pos[i]+1):(equal.pos[i+1]-2)]
    params[[param.name]] <- param.data
  }
  ## last param
  param.name <- raw.input[equal.pos[n.params]-1]
  param.data <- raw.input[(equal.pos[n.params]+1):length(raw.input)]
  params[[param.name]] <- param.data
  return(params)
}

read.chem.table <- function( table.name ) {
  orig.table <- read.table(table.name, header=T, sep="\t", comment.char="!")
  named.table <- orig.table[,2:ncol(orig.table)]
  rownames(named.table) <- orig.table[,1]
  return(named.table)
}

init.atom.mass <- function() {
  atom.mass.vec <- c(12.000000, #C
                     1.007825,  #H
                     15.994915, #O
                     14.003074, #N
                     31.972072, #S
                     30.973763, #P
                     15.000109, #N15
                     2.014102,  #H2
                     13.003355, #C13
					 17.999161, #O18
                     1.0072765,  #Hplus
                     34.96885, #Chlorine
                     78.91833, #Bromine
		     73.9224764 #Selenium
                     )
  names(atom.mass.vec) <- c("C","H","O","N","S","P","N15","H2","C13","O18","Hplus","Cl","Br","Se")
  return(atom.mass.vec)
}

init.aa.mass <- function(atom.mass.vec, chem.table ) {
  aa.names <- rownames(chem.table)
  chem.names <- colnames(chem.table)
  aa.mass.vec <- rep(0, length(aa.names))
  names(aa.mass.vec) <- aa.names
  for (aa in aa.names ) {
    mass <- 0.0
    for (chem in chem.names) {
      mass <- mass + atom.mass.vec[chem]*chem.table[aa,chem]
    }
    aa.mass.vec[aa] <- mass
  }
  return(aa.mass.vec)
}

element.count <- function(sequence.vec, element, chem.table) {
  ## sequence.vec is a vector of sequence character
  defined.aas <- rownames(chem.table)
  if ( element %in% colnames(chem.table) ) {
    count<- 0
    for ( aa in sequence.vec ) {
      if ( aa %in% defined.aas ){
        count<- count + chem.table[aa,element]
      }
    }
    count <- count+ chem.table["NTERM",element] + chem.table["CTERM",element]
    return(count)
  } else {
    return(0)
  }
}

vectorize.sequence <- function(sequence) {
  ## get rid of flanking residues deliminated by "."
  peptide.vec <- unlist( strsplit(sequence,".",fixed=T) )
  peptide.vec <- unlist( strsplit(peptide.vec[2],"",fixed=T) )
  return(peptide.vec)
}

calc.peptide.mass <- function(sequence, aa.mass.vec) {
  peptide.vec <- vectorize.sequence(sequence)
  mass <- 0
  defined.aas <- names(aa.mass.vec)
  for ( aa in peptide.vec ) {
    if ( aa %in% defined.aas ) {
      mass <- mass + aa.mass.vec[aa]
    }
  }
  mass <- mass + aa.mass.vec["NTERM"] + aa.mass.vec["CTERM"]
  return(mass)
}


calc.num.elements <- function( sequence, chem.table) {
  peptide.vec <- vectorize.sequence(sequence)
  elements <- colnames(chem.table)
  elements.count <- rep(0, length(elements) )
  names(elements.count) <- elements
  for ( e in elements) {
    elements.count[e] <- element.count(peptide.vec, e, chem.table)
  }
  return(elements.count)
}

calc.num.N15 <- function(sequence, chem.table) {
  peptide.vec <- vectorize.sequence(sequence)
  num.N15 <- element.count(peptide.vec, "N15", chem.table)
  return(num.N15)
}


