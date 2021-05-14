
#default.postcode.file = "Data/Postcode_Estimates_Table_with_coordinates2.csv"
source("utils.R")
library(sf)
library(KernSmooth)
library(spatstat)
library(truncnorm)
source("surveillance_utils.R")
source("plot_utils.R")

#load('Data/postcode2coord.Rdata')
# UK range (latitude and longitude)
#Y.range = range(postcode2coord[!is.na(postcode2coord$latitude),]$latitude) + c(-2,2) * var(postcode2coord[!is.na(postcode2coord$latitude),]$latitude)
# [1] 50.83559 55.64878
#X.range = range(postcode2coord[!is.na(postcode2coord$longitude),]$longitude) +  c(-2,2) * var(postcode2coord[!is.na(postcode2coord$longitude),]$longitude)
# [1] -4.717408  0.348140

ranScanCreateCylinders<-function(observation.matrix, baseline.matrix, emmtype,
                                 week.range, n.cylinders=1000, rs=0.1,
                                 p.val.threshold=0.05, plot=F, coord.df=postcode2coord,  size_factor=1){
  # n.cylinders=10000 takes around 3 hours for the whole dataset, to end up with 300 non-empty cylinders
  # observation.matrix  and baseline.matrix have dimension 
  # This function scans the matrices
  # -1 in observation.matrix index means that we are excluding from week NA
  #' @param emmtype
  #' @param week.range
  #' @param n.cylinders (integer) number of proposed cylinders per week interval.
  #' @param rs
  #' @param p.val.threshold
  #' @param plot (boolean)
  observation.matrix = observation.matrix[!(rownames(observation.matrix) == 'NA'),]
  baseline.matrix = baseline.matrix[!(rownames(baseline.matrix) == 'NA'),]
  week.range = range(as.integer(week.range))
  
  print(sum(baseline.matrix) )
  print(sum(observation.matrix))
  
  if (sum(baseline.matrix) > 2 * sum(observation.matrix)){
    print(c(sum(baseline.matrix),  sum(observation.matrix) ))
    writeLines("Warning: the baseline is too high. Have you multiplied for the emmtype factor?")
    writeLines("E.g., `ranScanCreateCylinders(observation.matrix, baseline.matrix*emmtype.factor[['12.0']], '12.0', starting.week, n.cylinders=100, rs=0.1)`")
  }
  init=Sys.time()
  coord.df = coord.df[!is.na(coord.df$latitude),]
  coord.km.df = coord.df
  coord.km.df[,2:3] = vlatlong2km(coord.df[,2:3])
  if ((week.range[1] < min(as.integer(colnames(observation.matrix)))) | (week.range[2] > max(as.integer(colnames(observation.matrix))))){
    # line = sprintf("Try with `starting.week` < %d", ncol(observation.matrix)-2)
    # writeLines(c("`starting.week` is bigger than the matrix length.", line))
    A = sprintf("%d-%d", week.range[1], week.range[2])
    B  = range(as.integer(colnames(observation.matrix)))
    B = sprintf("%d-%d", B[1], B[2])
    writeLines(paste0("Error: `week.range`` is ", A, ", while the range of `observation.matrix` is ", B, "." ))
    return(NA)
  }
  # weeks = as.integer(colnames(observation.matrix)[-seq(1:(starting.week+2))])
  # weeks = weeks[seq(1, length(weeks), 20)]
  # cylinders0 = data.frame(x=double(), y=double(), rho=double(), t.low=integer(), t.upp=integer(), n_obs=integer(), mu=double(), lower=integer(), upper=integer(), p.val=double(), warning=logical())
  radia_and_heights = f_radia_and_heights(baseline.matrix, 1:24) * size_factor
  cat("Evaluating exceedances from ",
        as.character(week2Date(week.range[1])),   # this interactive output doesnt work in the prospective mode
        " to ",
        as.character(week2Date(week.range[2])),        
        " for emm type ", emmtype, ".\n")
  
  # generate cylinders
  cylinders = rcylinder2(n.cylinders, observation.matrix, week.range, radia_and_heights, coord.km.df)
  if (NROW(cylinders) > 0){
    cylinders[,c('n_obs', 'mu', 'lower', 'upper', 'p.val')] = t(apply(cylinders, 1, compute,
                                                    observation.matrix, baseline.matrix, coord.km.df))
    ## da vettorizzare:
    # cylinders$warning = apply(cylinders, 1, function(x){ifelse((x['p.val'] < p.val.threshold) & (x['n_obs'] > 0), TRUE, FALSE)})
    # vettorizzato:
    cylinders$warning = cylinders['p.val'] < p.val.threshold
  }else{
    cat("No cases in the selected week range. No cylinder list returned.\n")
  }
  print(Sys.time() - init)
  return(cylinders)
}

# 
# if (plot){
#   output.basename = paste0("Exceed_", emmtype, "_", as.character(week2Date(week.range[1])), "_", as.character(week2Date(week.range[2])))
#   # For each location, count how many cylinders are warning and how many are not.
#   warning.ratio = sapply(1:nrow(coord.df),
#                          FUN=warning_ratio,
#                          observation.matrix[,week.range[1]:week.range[2]], cylinders, coord.df)
#   
#   idx = which(warning.ratio > p.val.threshold)
#   png(paste0(output.basename, '.png'), width = 3.25, height = 3.25, units = 'in', res=1200, pointsize = 4)
#   par(mar=c(2,2,4,2))
#   plot_map( X.range+c(1,-1), Y.range, main=paste("Exceedances from", week2Date(week.range[1]), "to", week2Date(week.range[2])))
#   plot_cylinders(cylinders[cylinders$warning,], lwd=0.4)
#   plot_cases(observation.matrix, week.range,  coord.df)
#   plot_cases(observation.matrix[idx,], week.range,  coord.df[idx,], col='red')
#   dev.off()
#   # save cases on file
#   tmp.table = observation.matrix[idx, week.range[1]:week.range[2]]
#   colnames(tmp.table) = sapply(colnames(tmp.table), function(x){as.character(week2Date(as.integer(x)))})
#   write.table(tmp.table,
#               file=paste0(output.basename, ".csv"), sep=',', row.names = T, col.names = T, quote = T)
# }

ranScanPlotCylindersCI<-function(cylinders, confidence.level=0.68, title=NULL){
  plot(n_obs ~ mu, cylinders, title=title, xlab = 'Expected cases', ylab = 'Observed cases', pch=20, cex=0.8)
  mu=max(cylinders$mu)
  lines(0:mu, 0:mu)
  x = seq(0,mu,0.1)
  upper=qpois(confidence.level, x)
  lower=qpois(1-confidence.level, x)
  polygon(c(x, x[length(x):1]),c(upper, lower[length(lower):1]), col="#a6a6a666", border = "#a6a6a600")
  points(cylinders[cylinders$warning,]$mu, cylinders[cylinders$warning,]$n_obs, col='red', pch=20, cex=0.9)
}

ranScanSaveCylinders<-function(cylinders, file.basename){
  write.csv(cylinders, paste0(file.basename,'.csv'), quote = F, row.names = F)
  save.and.tell("cylinders", file = paste0(file.basename,'.Rdata'))
}

ranScanEvaluate<-function(case.file, cylinders, emmtype, p.val.threshold=0.05,
                          warning.score.name='warning.score'){
  case.df<-tryCatch({
    load(paste0(case.file, ".Rdata"))
    case.df
  },
  error = function(e){
    case.df = ranScanInit(case.file)
    return(case.df$case.df)
  })
  # if (!is.null(p.val.threshold)){
  #   cylinders$warning = apply(cylinders, 1, function(x){ifelse((x['p.val'] < p.val.threshold) & (x['n_obs'] > 0), TRUE, FALSE)})
  # }
  if (!('x' %in% names(case.df) & ('y' %in% names(case.df)))){
    writeLines("Inserting coordinates...")
    if (!('latitude' %in% names(case.df) & ('longitude' %in% names(case.df)))){
      case.df[,c("latitude", "longitude")] = t(apply(case.df, 1, postcode.to.location2))      
    }
    case.df[,c("y", "x")] = vlatlong2km(case.df[,c("latitude", "longitude")])
    writeLines("...Done.")
    save.and.tell('case.df', file=paste0(case.file, ".Rdata"))
  }
  idx = case.df$emmtype == emmtype
  writeLines(paste0("Computing warning scores for emmtype ", emmtype, "..."))
  case.df[idx,warning.score.name] = apply(case.df[idx,], 1, FUN=warning.score, cylinders)
  writeLines("...Done.")
  return(case.df)  
}


PostcodeAreasBoundaries=st_read("shape_file_new/Areas.shp")
# Keep only Englandpostcode areas:
PostCodeAreas=as.character(PostcodeAreasBoundaries[[1]])
PostCodeAreasC = read.csv("Data/Book1.csv")
PostCodeAreasC = PostCodeAreasC[PostCodeAreasC$Country=='England', ]
id = sapply(PostCodeAreas, function(x){x %in% PostCodeAreasC[,1]})
PostcodeAreasBoundaries = st_geometry(PostcodeAreasBoundaries)[id]
rm(id)
RegionBoundaries = st_read("shape_files/RegionBoundaries.shp")
RegionBoundaries = st_transform(RegionBoundaries, crs=st_crs(PostcodeAreasBoundaries))
UK = st_read("Other maps/Map_UK.shp")
UK = st_transform(UK, crs=st_crs(PostcodeAreasBoundaries))
Wales. = st_union(UK[UK$NAME_1=='Wales',])

plotBaseMap<-function(main=NULL, add=T, Wales=TRUE, ...){
  plot(st_geometry(UK), col=NA, border='#e2e2e288', lwd=0.8,  add=add,...)
  if (Wales){
    plot(st_geometry(Wales.), axes=F, border='black', lwd=1, add=T)    
  }
  plot(st_geometry(RegionBoundaries), axes=F, border='black', lwd=1,  main=main, add=T)
  # box(lwd=0.1)
}
X.range = c(-6,2)
Y.range = c(50,56)
#EnglandBoundaries=st_union(st_geometry(RegionBoundaries))

is_inside<-function(i,j,x,y){
  latitude=y[i]
  longitude=x[j]
  stp=st_point(c(longitude, latitude))
  return(!is.empty(st_contains(st_geometry(EnglandBoundaries), stp)[[1]]))
}


plotBaseline<-function(week, emmtype, z, add=F, tf=factor, emmf=emmtype.factor){
  xmax=max(z$fhat) * tf[week] * emmf[emmtype][[1]]
  ColorRamp=colorRampPalette(c("white", blues9))(100)
  ColorRamp_ex=ColorRamp[1:(xmax / 3 * 100)]
  image(z$x1, z$x2,
        z$fhat * tf[week] *emmf[emmtype][[1]],
        col=ColorRamp_ex, add = add, axes=F, xlab=NA, ylab=NA)
}

# ranScanMovie<-function(cylinders, observation.matrix, postcode2coord, output.basename, emmtype){
#   idx = rownames(observation.matrix) != 'NA'
#   observation.matrix = observation.matrix[idx,]
#   idx = rownames(postcode2coord) != 'NA'
#   postcode2coord = postcode2coord[idx,]
#   z = bkde2D(cbind(postcode2coord$longitude, postcode2coord$latitude), bandwidth = 0.1, gridsize = c(400,400))
#   if (dim(postcode2coord)[1] != dim(observation.matrix)[1]){
#     writeLines("`postcode2coord` and observation matrices must have the same size")
#     return(0)
#   }
#   
#   cylinders[,c(2,1)] = vkm2latlong(cylinders[,c(1,2)])
#   
#   cRamp= colorRamp(c('Blue','Red'))
#   for(i in 2:ncol(observation.matrix)){
#     cat("week",i,"\r")
#     A=paste0(output.basename, '_', i, '_' , emmtype, '_.png')
#     png(A, width = 3.25, height = 3.25, units = 'in', res=600, pointsize = 4)
#     par(mar=c(2,2,4,2))
#     main=paste0("week ", i-2, ", emmtype ", emmtype)
#     plotBaseline(i, emmtype, z, add=F)
#     plotBaseMap(main='', add=T)
#     text(x = -3.5, y = 52.5, main) 
#     idx = observation.matrix[,i] > 0
#     warning.color = rgb2hex(cRamp(
#       sapply(1:nrow(postcode2coord[idx,]),
#              FUN=warning_ratio2, observation.matrix[idx,i],i, cylinders, postcode2coord[idx,])))
#     points(postcode2coord[idx,]$longitude, postcode2coord[idx,]$latitude, pch=20, cex=1.5, col=warning.color)
#     title(week2Date(i))
#     dev.off()
#   }
# #  shell("magick Fig/*.png Fig/movie.gif")
# }


# ranScanPlotCluster<-function(observation.matrix, postcode2coord.km, main=''){
#   library(tsne)
#   idx = rownames(observation.matrix) != 'NA'
#   observation.matrix = observation.matrix[idx,]
#   idx = rownames(postcode2coord.km) != 'NA'
#   postcode2coord.km = postcode2coord.km[idx,]
#   cases = which(observation.matrix>0, arr.ind=TRUE)
#   cases = as.data.frame(cases)
#   c3 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 3]})
#   c2 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 2]})
#   cases2 = cbind(cases,c2,c3) 
#   cases2 = cases2[!is.na(cases2[,3]), ]
#   X = tsne(cases2[,2:4])
#   palette = colorRampPalette(c('blue', 'red'))(max(cases2[,2]))
#   warning.marker = ifelse(8, 20, sapply(1:nrow(postcode2coord[idx,]),
#                                         FUN=warning_ratio2, observation.matrix[idx,i],i, cylinders, postcode2coord[idx,])>0.9 )
#   
#   plot(X[,1],X[,2], xlab='C1', ylab='C2', col=colormap[cases2[,2]], pch=warning.marker, main=main)
# }


ranScanCluster<-function(case.df, emmtype, makeplot=FALSE, warning.score='warning.score'){
  library(tsne)
  # idx = rownames(observation.matrix) != 'NA'
  # observation.matrix = observation.matrix[idx,]
  # idx = rownames(postcode2coord.km) != 'NA'
  # postcode2coord.km = postcode2coord.km[idx,]
  # cases = which(observation.matrix>0, arr.ind=TRUE)
  # cases = as.data.frame(cases)
  # c3 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 3]})
  # c2 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 2]})
  # cases2 = cbind(cases,c2,c3)
  idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
  case.df = case.df[idx, ]
  X = tsne(case.df[,c('x', 'y', 'SAMPLE_DT_numeric')])
  # palette = colorRampPalette(c('blue', 'red'))(max(case.df$SAMPLE_DT_numeric))
  if (makeplot == TRUE){
    library(viridis)
    palette = viridis(max(case.df$SAMPLE_DT_numeric), begin = 0, end = 1)
    warning.marker = 20 #ifelse(case.df[,warning.score] > 0.9, 20, 1)
    plot(X[,1],X[,2], xlab='C1', ylab='C2', col=palette[case.df$SAMPLE_DT_numeric], pch=warning.marker,
         main=emmtype)
    
    palette = viridis(max(case.df$SAMPLE_DT_numeric), alpha = 1, begin = 0, end = 1)
    idx = ifelse(case.df[,warning.score] > 0.9, T, F)
    points(X[idx,1],X[idx,2], xlab='C1', ylab='C2', col=palette[case.df[idx, "SAMPLE_DT_numeric"]], pch=1, cex=2)
  }
  return(X)
}

ranScanPlotCluster0<-function(X, case.df, emmtype, warning.score='warning.score', ...){
  threshold = 0.95
  idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
  case.df = case.df[idx, ]
  palette = colorRampPalette(c('blue', 'red'))(101)
  par(fig=c(0, 1, 0, 1))
  par(mar=c(4, 4, 0.9, 4))
  #main=paste0('Embedding for emmtype ',emmtype)
  plot(X[,1], X[,2], xlab='C1', ylab='C2', col=palette[round(case.df[,warning.score] * 100+1)], pch=19, ...)
  idx2 = (case.df[,warning.score]>0.95)
  points(X[idx2,1], X[idx2,2], cex = 2)
  
  legend("bottomleft", c("cases", as.expression(bquote(italic(w) ~ ">" ~ .(threshold)))),
         col = c(palette[as.integer(length(palette)/2)], 'black'), pch=c(19,1), pt.cex = c(1, 2))
  
  
  Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
  par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
  par(mar=c(0,0,0,0))
  plot(c(0, 1), c(0, 1), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
  rasterImage(Cb[seq(length(palette), 1, -1),], 0, 0, 1, 1)
  rect(0, 0, 1, 1)
  axis(side = 4)
  mtext("warning score", side = 4, line = 2)
}

ranScanPlotCluster<-function(X, case.df, emmtype, warning.score='warning.score', threshold=NULL, ...){
  if (is.null(threshold)){
    ranScanPlotCluster0(X, case.df, emmtype, warning.score, ...)
  }else{
    library(viridis)
    idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
    case.df = case.df[idx, ]
    palette = viridis(max(case.df$SAMPLE_DT_numeric)+1, alpha = 1, begin = 0, end = 1)
    par(fig=c(0, 1, 0, 1))
    par(mar=c(4, 4, 0.9, 4))
    plot(X[,1],X[,2], xlab='C1', ylab='C2', col=palette[case.df$SAMPLE_DT_numeric+1], pch=19, ...)
    # mtext(paste0('Embedding for emmtype ',emmtype), 3)
    palette = viridis(max(case.df$SAMPLE_DT_numeric), alpha = 1, begin = 0, end = 1)
    idx = ifelse(case.df[,warning.score] > threshold, T, F)
    points(X[idx,1], X[idx,2], xlab='C1', ylab='C2', col='black', pch=1, cex=2)
    legend("bottomleft", c("cases", as.expression(bquote(italic(w)~ ">" ~ .(threshold)))),
           col = c(palette[as.integer(length(palette)/2)], 'black'), pch=c(19,1), pt.cex = c(1, 2))
    
    Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
    par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
    par(mar=c(0,0,0,0))
    plot(c(0, 1), c(1, len(palette)), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
    rasterImage(Cb[seq(length(palette), 1, -1),], 0, 0, 1, len(palette))
    rect(0, 0, 1, length(palette))
    rr = max(case.df$SAMPLE_DT_numeric)
#    ticks = as.integer(seq(1, rr, (rr - 1)/ 10))
    ticks = c(0, 53, 53*2-1, 53*3-2, 53*4-3, 53*5-4)
    axis(side = 4, at=ticks, labels=week2Year(ticks))
#    mtext("time [week]", side = 4, line = 2)
  }
}

ranScanPlotCluster3<-function(X, case.df, emmtype,
                              warning.score='warning.score',
                              threshold=0.95, Legend = T,...){
    idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
    case.df = case.df[idx, ]
    ll = c("East Midlands", "East of England", "London", "North East",
           "North West", "South East", "South West", "West Midlands","Yorkshire&Humber")
    
    if(is.null(case.df$Postcode.Region)){
      load("Data/Area2Region_list.RData")
      
      case.df$Postcode.Region = ordered(unlist(sapply(case.df$`Patient Postcode`,
                                                     postcode.to.region,
                                                     Area2Region_list)), levels = ll)
    }
    palette = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#0000FF")
    palette.alpha = c("#1B9E7755", "#D95F0255", "#7570B355",
                      "#E7298A55", "#66A61E55", "#E6AB0255",
                      "#A6761D55", "#66666655", "#0000FF55")
    
    case.df$color.Region = palette[case.df$Postcode.Region]
    case.df$color.Region.alpha = palette.alpha[case.df$Postcode.Region]
    
    par(fig=c(0, 1, 0, 1))
    par(mar=c(4, 4, 0.9, 4))
    
#    paste0('Embedding for emmtype ',emmtype)
#    palette = viridis(max(case.df$SAMPLE_DT_numeric), alpha = 0.1, begin = 0, end = 1)
    plot(X[,1],X[,2], xlab='C1', ylab='C2',
         col=case.df$color.Region.alpha, pch=19, ...)


    idx2 = ifelse(case.df[,warning.score] > threshold, T, F)
    while((sum(idx2) == 0) & (threshold > 0)){
      threshold = threshold - 0.1
      print(threshold)
      idx2 = ifelse(case.df[,warning.score] > threshold, T, F)
    }
    
    region = ordered(case.df[idx2,]$Postcode.Region)
    as.int = as.integer(region)

    points(X[idx2,1], X[idx2,2], xlab='C1', ylab='C2', col= case.df[idx2,]$color.Region,
           pch=1, cex=1.8)
    
    # if ((length(levels(region)) > 0) & (Legend == T)){
    #   legend("bottomleft", levels(region),
    #          col = palette[1:length(levels(region))],
    #          pch = 1, cex=0.7, pt.cex = 1.4)
    # }
    
    legend('bottomleft',
           c(as.character(sort(unique(case.df$Postcode.Region))),
             as.expression(bquote(italic(w) ~ '>' ~ .(threshold)))),
           col=c(palette.alpha[sort(unique(case.df$Postcode.Region))], 'black'),
           pch=c(rep(19, len(palette[sort(unique(case.df$Postcode.Region))])), 1),
           pt.cex = c(rep(1, len(palette[sort(unique(case.df$Postcode.Region))])), 1.4),
           cex=0.7
           )
}

#' ranScanMutateCylinders<-function(observation.matrix, baseline.matrix, cylinders, n.cylinders=NULL,
#'                                  plot=F, p.val.threshold=0.05){
#'   # observation.matrix  and baseline.matrix have dimension 
#'   # This function scans the matrices 
#'   # -1 in observation.matrix index means that we are excluding from week NA
#'   #' @param emmtype
#'   #' @param starting.week
#'   #' @param n.cylinders
#'   #' @param rs
#'   #' @param p.val.threshold
#'   #' @param plot (boolean)
#'   #' 
#'   if (sum(baseline.matrix) > 2 * sum(observation.matrix)){
#'     writeLines("Warning: the baseline is too high. Have you multiplied for the emmtype factor?")
#'   }
#'   init=Sys.time()
#'   if (is.null(n.cylinders)){
#'     n = nrow(cylinders)    
#'   }else{
#'     n=n.cylinders
#'   }
#'   idx = sample(1:nrow(cylinders), size=n, prob = cylinders$n_obs, replace = T)
#'   cylinders = cylinders[idx,]
#'   cylinders$x = cylinders$x + rnorm(n, sd=0.3)
#'   cylinders$y = cylinders$y + rnorm(n, sd=0.3)
#'   b = round(max(cylinders$t.upp) - cylinders$t.upp)
#'   shift = floor(rtruncnorm(n, b=b, mean=0, sd=2))
#'   cylinders$t.low = cylinders$t.low + shift
#'   cylinders$t.upp = cylinders$t.upp + shift
#'   
#'   cylinders[,c('n_obs', 'mu', 'lower', 'upper', 'p.val')] = t(apply(cylinders, 1, compute,
#'                                                                     observation.matrix, 
#'                                                                     baseline.matrix,
#'                                                                     postcode2coord))
#'   cylinders$warning = apply(cylinders, 1, function(x){ifelse((x['p.val'] < p.val.threshold) & (x['n_obs'] > 0), TRUE, FALSE)})
#'   print(Sys.time() - init)
#'   return(cylinders)
#' }

