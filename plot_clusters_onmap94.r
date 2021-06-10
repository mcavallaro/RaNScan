
emmtype = '94.0'
# source("ranScanInit.r")
case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"

ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'SAMPLE_DT_numeric')
#ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'RECEPT_DT_numeric')

load(paste0("~/Outbreak/Data/", emmtype, "_obs.Rdata"))
baseline.matrix = ranScanCreateBaselineMatrix(case.file)
emmtype.factor = ranScanEmmtypeFactor(case.file)

n_cyl = 20000
week.range = c(1, 298)

set.seed(1)
cylinders.94 = ranScanCreateCylinders(observation.matrix,
                                   baseline.matrix * emmtype.factor[[emmtype]], emmtype, week.range,
                                   n_cyl, postcode2coord,  size_factor = 1.1)
set.seed(1)
case.df.tmp.94 = ranScanEvaluate(case.file, cylinders.94, emmtype, p.val.threshold=0.05)

idx = ((case.df.tmp.94$emmtype == emmtype) & !is.na(case.df.tmp.94$x))
case.df.tmp1 = case.df.tmp.94[idx, ]
set.seed(1); tmp.94 = ranScanCluster(case.df.tmp1, emmtype)


warning.threshold = 0.95

png(paste('Manuscript/', emmtype, '.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.94, case.df.tmp1, emmtype)
dev.off()

png(paste('Manuscript/', emmtype, 'thresh.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.94, case.df.tmp1, emmtype, threshold=warning.threshold, xlim=c(-21, 18))
dev.off()

png(paste('Manuscript/', emmtype, 'regions.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster3(tmp.94, case.df.tmp1, emmtype, Legend = T, xlim=c(-32, 18), threshold = 0.95)
# legend('bottomleft', c('North West, w>0.85','South East, w>0.85','South West, w>0.85'),
#        col=c("#1B9E77", "#D95F02", "#7570B3"), pch=1, pt.cex=c(1.4, 1.4, 1.4), cex=0.7)
dev.off()


# idx.a = tmp.94[,2] > 10
# idx.b = case.df.tmp1$warning.score > 0.75
# idx.c = idx.a & idx.b
# 
idxw = case.df.tmp1$warning.score > 0.95
idxw = idxw & tmp.94[,1] > 7 & tmp.94[,2] > -32

palette = plasma(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)+1  - min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric), begin = 0, end = 0.75)

png('Manuscript/_94.0clust.png', res=500, width = 5, height = 4, units='in')

par(fig=c(0, 1, 0, 1))
par(mar=c(0, 0, 0, 4))

case.df.tmp1[idxw,]$latitude %>% range
case.df.tmp1[idxw,]$longitude %>% range

plotBaseMap(add=F, xlim=c(-6,0), ylim=c(50.5, 51)) #emm type 94

colors = palette[case.df.tmp1[idxw,]$SAMPLE_DT_numeric - min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric) + 1]
points(latitude ~ longitude,
       data=case.df.tmp1[idxw,],
       col=colors, cex=1.8)
#
# points(latitude ~ longitude,
#        data=case.df.tmp1[idx.d & idxw,],
#        col='red', cex=1.8, pch=4)

# mtext(paste0("emmtype ", emmtype), side = 1)

par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
dev.off()
