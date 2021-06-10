
# INITIALISE:
source("ranScanInit.r")

case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"
# 
# tmp = ranScanInit(case.file)
# ranScanCreateObservationMatrices(case.file, tmp$emmtypes[1:4])
# 
# postcode2coord = ranScanPostcodeMap(observation.matrix)
# load("Data/1.0_obs.Rdata")
# baseline.matrix = ranScanCreateBaselineMatrix(case.file)
# emmtype.factor = ranScanEmmtypeFactor(case.file)
# time.factor = ranScanTimeFactor(case.file)
# 
# #
# sum(baseline.matrix * emmtype.factor[['1.0']], na.rm=T)
# sum(observation.matrix, na.rm=T)
# 
# load("Data/12.0_obs.Rdata")
# sum(baseline.matrix * emmtype.factor[['12.0']], na.rm=T)
# sum(observation.matrix, na.rm=T)
# 
# source("ranScanEvaluate.r")
# ranScanCreateObservationMatrices(case.file, c('33.0', '108.1', '12.0'))
# ranScanCreateObservationMatrices(case.file, c('44.0'))
#for (n_cyl in c(20000, 30000)){

A = list()
n_cyl = 30000
for(sizefactor in c(1)){
  for (emmtype in c('94.0', '44.0', '12.0', '33.0', '108.1')){
    t1 = Sys.time()
    load(paste0("~/Outbreak/Data/", emmtype, "_obs.Rdata"))
    cylinders = ranScanCreateCylinders(observation.matrix,
                                       baseline.matrix * emmtype.factor[[emmtype]], emmtype, week.range,
                                       n_cyl, postcode2coord,  size_factor = sizefactor)
    
    
    case.df.tmp = ranScanEvaluate(case.file, cylinders, emmtype, p.val.threshold=0.05)
    
    A[[emmtype]] = case.df.tmp[case.df.tmp$emmtype == emmtype, ]
    tmp = ranScanCluster(case.df.tmp, emmtype)
    
    png(paste('Manuscript/n', as.character(n_cyl),as.character(sizefactor), emmtype, '.png', sep='_'),
        res=600, width=5, height=4, units="in")
    ranScanPlotCluster(tmp, case.df.tmp, emmtype)
    dev.off()
    
    png(paste('Manuscript/n', as.character(n_cyl),as.character(sizefactor), emmtype, 'thresh.png', sep='_'),
        res=600, width=5, height=4, units="in")
    ranScanPlotCluster(tmp, case.df.tmp, emmtype, threshold=0.9)
    dev.off()
    
    png(paste('Manuscript/n', as.character(n_cyl), as.character(sizefactor), emmtype, 'regions.png', sep='_'),
        res=600, width=5, height=4, units="in")
    ranScanPlotCluster3(tmp, case.df.tmp, emmtype, threshold=0.7, Area2region_list=Area2Region_list)
    dev.off()
    
    
    png(paste('Manuscript/n', as.character(n_cyl), as.character(sizefactor), emmtype, 'clustmap.png', sep='_'),
        res=600, width=5, height=4, units="in")
    plotBaseMap(main='', add=T)
    points
    dev.off()
    # ranScanSaveCylinders(cylinders, paste0('Data/cylinders_for_', emmtype, '.csv'))
    # mutcylinders = ranScanMutateCylinders(observation.matrix, baseline.matrix* emmtype.factor[[emmtype]], cylinders, n.cylinders=40000)
    # ranScanSaveCylinders(mutcylinders, paste0('Data/mutcylinders_for_', emmtype, '.csv'))
    # case.df = ranScanEvaluate(case.file, cylinders, emmtype, p.val.threshold=0.05)
    # write.csv(case.df[case.df$emmtype==emmtype,], file = paste0('Data/cases_', emmtype, '_with_warnings.csv'))
    # write.csv(case.df,file = paste0('Data/cases_full_', emmtype, '_with_warnings.csv'))
    cat("n_cyl", n_cyl, "emmtype", emmtype)
    print(Sys.time() - t1)
  }
}
# 
# 
# emmtype = '33.0'
# load(paste0("~/Outbreak/Data/", emmtype, "_obs.Rdata"))
# load(paste0("~/Outbreak/Data/mutcylinders_for_",emmtype,".csv.Rdata"))
# cidx = !(row.names(observation.matrix) == 'NA')
# ranScanMovie(cylinders, observation.matrix[idx,], postcode2coord, 'Fig/2Dec', emmtype)



# East of England 33	         ? - 18-08-2019  44
# South West 21	30-01-2017 - 05-03-2018  94
# Midlands 18	26-11-2018 - 07-04-2020 108.1
# West Yorkshire 36	23-02-2018 - 07-05-2019 108.1
# 
# date2week<-function(date){
#   date = as.POSIXct(paste0(as.character(date), " UTC"), tz = "UTC")
#   return(
#     as.integer(difftime(date, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
#   )
# }
# 
# 
# cat('-- ', date2week("2019-08-18"))
# cat(date2week("2017-01-30"), date2week("2018-03-05"))
# cat(date2week("2018-11-26"), date2week("2020-04-07"))
# cat(date2week("2018-02-23"), date2week("2019-05-07"))
# # data,range emmtype
# # --  241 44
# # 108-165 94
# # 203-274 108.1
# # 164-226 108.1
# 
png('Manuscript/prevalence.png', res=500, width=5 * 2, height=4, units="in", pointsize = 12.1)
par(mfrow=c(1, 2), mar=c(5,4.3,4,2) + 0.1)
plot.emm.fraction(case.df , c('1.0', '89.0', '12.0'), 'topleft', ylim=c(0,0.51), xlab=expression(tau), cex.lab=1.5)
mtext('A', 3, at = max(case.df$SAMPLE_DT), cex=2, padj=1.5)
plot.emm.fraction(case.df , c('94.0','108.1','44.0', '33.0'), 'topleft', ylim=c(0,0.06), xlab=expression(tau), cex.lab=1.5)
mtext('B', 3, at = max(case.df$SAMPLE_DT), cex=2, padj=1.5)
dev.off()



png('Manuscript/time_factor2.png', res=600, width=5, height=4, units="in", pointsize = 10)
time.series = table(case.df$SAMPLE_DT_numeric)
plot(time.series[-1], xlab = 't', ylab='No. of cases / week', xaxt='n', type = 's')
axis(1,  at=0:300, NA, tck = -0.035, lwd=0.01) 
axis(1, at=c(0,53, 105,  157,209,261),  week2Year(c(0,53, 105,  157,209,261)), tck = -0.07) 
x <- as.numeric(dimnames(observation.matrix)[[2]][-1])
#time.factor = ranScanTimeFactor(case.file, F)
ci = predict.cmle.ci(x, attributes(time.factor.tmp1)$Parameters)
col = paste0(col2hex('red'), '22')
polygon(c(x, x[length(x):1]), c(ci$upper, ci$lower[length(ci$lower):1]), border=NA, col=col)
lines(x, predict.cmle(x, attributes(time.factor.tmp1)$Parameters), col='red', lwd=2, lty=5)
legend("topleft", pch=c(NA,NA), lty = c(1, 5), lwd=c(1, 1.5), fill=c(NA, col), border = c(NA, NA),
       col=c('black', 'red'), c("Observation", "Baseline"))
dev.off()
