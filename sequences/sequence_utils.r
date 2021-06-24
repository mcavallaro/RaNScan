#fileName <- 'foo.txt'
#readChar(fileName, file.info(fileName)$size)    

dna_string="GATTACAGATTACA"


dna2bitstring<-function(dna_string){
  splitted = strsplit(dna_string, '')[[1]]
  bitstring = vector(length = length(splitted) * 2, mode = 'logical')
  for (i in 1:length(splitted)){
    b = c(i,(i+1)) + (i-1)
    bitstring[b] = switch(splitted[i], 'G'=c(F,T), 'A'=c(F,F), 'T'=c(T,T), 'C'=c(T,F),
                                       'g'=c(F,T), 'a'=c(F,F), 't'=c(T,T), 'c'=c(T,F))
  }
  return(bitstring)
}

distance_ape<-function(a,b,...){
  m = matrix(c(strsplit(a,'')[[1]], strsplit(b,'')[[1]]), nrow = 2, byrow = T)
  return(dist.gene(m,...)[1])
}


hamming_distance0<-function(bitstring1, bitstring2){
  pseudoH = bitstring1 == bitstring2
  L = length(bitstring1)
  pseudoHa = pseudoH[seq(1,L,2)]
  pseudoHb = pseudoH[seq(2,L,2)]
  H=L/2 - sum(pseudoHa & pseudoHb)
  return(H)
}



hamming_distance<-function(dna_string1, dna_string2){
  bitstring1=dna2bitstring(dna_string1)
  bitstring2=dna2bitstring(dna_string2)
  return(hamming_distance0(bitstring1, bitstring2))
}


