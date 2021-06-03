source("sequence_utils.R")
file.list = Sys.glob("hamming*.sds")
idx = order(global.scores$global.score, decreasing = T)
global.scores =  global.scores[idx,]
global.scores$emmtype = as.character(global.scores$emmtype )

part.file.list = c()
for (i in 1:length(file.list)){
  file = file.list[i]
  sub.string = substring(strsplit(file, '.sds', fixed=T)[[1]], 15)
  if (sub.string %in% global.scores$emmtype){
    part.file.list = c(part.file.list, file)
  }
}

A = list()
for (file.name in part.file.list){
  emmtype=substring(strsplit(file.name, '.sds', fixed=T)[[1]], 15)
  A[[emmtype]] = readChar(file.name, file.info(file.name)$size)  
}

global.scores$in.file = apply(global.scores, 1, function(x){if (x['emmtype'] %in% names(A)){T}else{F} })
  
H = matrix(0, ncol=sum(global.scores$in.file), nrow=sum(global.scores$in.file))
for (i in 1:ncol(H)){
  for (j in 1:nrow(H)){
    if (i < j){
      string.i = A[[global.scores[global.scores$in.file,][i, 'emmtype']]]
      string.j = A[[global.scores[global.scores$in.file,][j, 'emmtype']]]
      h = hamming_distance(string.i, string.j)
      H[i,j] = h
      H[j,i] = h
    }
  }
}

#par(mar=c(0,0,0,0))
#plot(c(1, 55), c(1, 55), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
image(1:51, 1:51, H / max(H), useRaster=F, axes=F)
axis(at = 1:51, side = 2)
axis(at = 1:51, side = 1)
