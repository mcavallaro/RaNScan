library(sf)
library(spatstat)


plot_cases<-function(observation.matrix, week.range, postcode.locations, col='#ffa500', cex=0.8, pch=NULL){
  weeks<-week.range[1]:week.range[2]

  case_postcodes = apply(observation.matrix[,weeks], 1 , function(x)any(as.logical(x)))
  
  cex = apply(observation.matrix[case_postcodes,weeks], 1 , function(x){sum(x)*cex})
  
  if (is.null(pch)){
    pch=16
  }
  case_locations = postcode.locations[case_postcodes,2:3]
  points(case_locations$longitude, case_locations$latitude, col=col, pch=pch, cex=cex)
}




plot_cylinders<-function(cylinders, color="#3182bd99",lwd=0.1){ #'#ffebeb88'
  for(i in 1:nrow(cylinders)){
    plot.owin(disc(cylinders$rho[i], c(cylinders$x[i],cylinders$y[i])), add=TRUE, col=NA, lwd=lwd, border=color)
  }
}


scale2one<-function(x){
  x = order(x)
  range = range(x)
  return( (x - range[1]) / (range[2] - range[1]))
  # x = (x - mean(reference)) / sqrt(sd(reference))  + 0.5
  # x = ifelse(x < 0, 0, x)
  # x = ifelse(x > 1, 1, x)
  # return(x)
}


plot.emm.fraction<-function(case.df, emmtypes, legend.position='topright', ...){
  case.df.ordered = case.df[order(case.df$SAMPLE_DT_numeric),]
  case.df.ordered$cum.cases = 1:nrow(case.df)
  ylab = expression(lambda[m])
  y = 0
  lwd = 2
  l = length(emmtypes)
  cols = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  for (i in 1:l){
    emmtype = emmtypes[i]
    y = c(y, range(cumsum(case.df.ordered$emmtype == emmtype)  / case.df.ordered$cum.cases))
  }
  x = case.df.ordered$SAMPLE_DT
  plot(range(x), range(y), type = 'n', ylab=ylab, ...)
  for (i in 1:l){
    emmtype = emmtypes[i]
    y = cumsum(case.df.ordered$emmtype == emmtype)  / case.df.ordered$cum.cases
    lines(x, y, col=cols[i], lwd=lwd)
  }
  legend(legend.position, lwd=lwd, col=cols[1:l], emmtypes)
}



plot.spaghetti<-function(prospective.w.df, spaghetti.col='black', ...){
  cols = 16:NCOL(prospective.w.df)
  x = 1:len(cols)
  plot(range(x), c(0.5,1), type='n', xaxt='n', ...)
  for (i in 1:2300){
    if( !is.na(prospective.w.df[i,'x'])){
      lines(x, unlist(prospective.w.df[i, cols]), col=spaghetti.col)
    }
  }
  axis(side = 1, at = x, labels = colnames(prospective.w.df)[cols])
  abline(h=0.95)
}


