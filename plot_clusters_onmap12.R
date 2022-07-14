source('ranScanInit2.R')
source('ranScanEvaluate.R')
source("plotBaseMap.r")
len=length

emmtype = '12.0'

case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"

ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'SAMPLE_DT_numeric')
#ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'RECEPT_DT_numeric')
load(paste0("Data/", emmtype, "_obs.Rdata"))

baseline.matrix = ranScanCreateBaselineMatrix(case.file)
emmtype.factor = ranScanEmmtypeFactor(case.file)

set.seed(1)
cylinders.12 = ranScanCreateCylinders(observation.matrix,
                                   baseline.matrix * emmtype.factor[[emmtype]], emmtype, week.range,
                                   n_cyl, postcode2coord,  size_factor = 1)
set.seed(1)
case.df.tmp.12 = ranScanEvaluate(case.file, cylinders.12, emmtype, p.val.threshold=0.05)

idx = ((case.df.tmp.12$emmtype == emmtype) & !is.na(case.df.tmp.12$x))
case.df.tmp1 = case.df.tmp.12[idx, ]

set.seed(1); tmp.12 = ranScanCluster(case.df.tmp1, emmtype)

png(paste('Manuscript/', emmtype, 'regions.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster3(tmp.12, case.df.tmp1, emmtype)
dev.off()


png(paste('Manuscript/', emmtype, '.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.12, case.df.tmp1, emmtype)
dev.off()

png(paste('Manuscript/', emmtype, 'thresh.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.12, case.df.tmp1, emmtype, threshold=0.95)
dev.off()

dev.off()
# idxw = case.df.tmp1$warning.score > 0.85
# plot(tmp.12[,1],tmp.12[,2])
# points(tmp.12[idxw,1],tmp.12[idxw,2],col='red')
# idxx = (tmp[,1] < 9.5) & (tmp[,1] > 3.5)
# idxy = (tmp[,2] < -3.5) & (tmp[,2] > -9.5)

# idxx1 = (tmp[,1] < 5) & (tmp[,1] > 3)
# idxy1 = (tmp[,2] < 3.5) & (tmp[,2] > 0)
# 
# idxx2 = (tmp[,1] < 10) & (tmp[,1] >4)
# idxy2 = (tmp[,2] < 10) & (tmp[,2] > 3.5)
# idxx = idxx1 | idxx2
# idxy = idxy1 | idxy2
# 
# idxw = case.df.tmp1[idxx & idxy,]$warning.score > 0.95
# idxw1 = case.df.tmp1[idxx1 & idxy1,]$warning.score > 0.95
# idxw2 = case.df.tmp1[idxx2 & idxy2,]$warning.score > 0.95
palette = plasma(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric) -
                   min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric), begin = 0, end = 0.75)

print(c(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric),
        min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)))

png('Manuscript/_12.0.clust.png', res=500, width = 5, height = 4, units='in')

par(fig=c(0, 1, 0, 1))
par(mar=c(0, 0, 0, 4))

case.df.tmp1[idxw,]$latitude %>% range
case.df.tmp1[idxw,]$longitude %>% range


plotBaseMap(add=F, xlim=c(-4.2,-0.4), ylim=c(50,54)) #emm type 94

colors = palette[case.df.tmp1[idxw,]$SAMPLE_DT_numeric - min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)]
points(latitude ~ longitude,
       data=case.df.tmp1[idxw,],
       col=colors, cex=1.7)
# colors = palette[case.df.tmp1[idxx2&idxy2,][idxw2,]$SAMPLE_DT_numeric - min(case.df.tmp1[idxx&idxy,][idxw,]$SAMPLE_DT_numeric)]
# points(latitude ~ longitude,
#        data=case.df.tmp1[idxx2 & idxy2,][idxw2,],
#        col=colors, cex=1.7)

#mtext(paste0("emmtype ", emmtype), side = 1)
par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
par(mar=c(0,0,0,0))
plot(c(0, 1), c(120, len(palette)+120), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
box(lwd=0.1)
Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
rasterImage(Cb[seq(length(palette), 1, -1),], 0, -0.04 * len(palette) + 120, 1, 1.035 * len(palette) + 120)
axis(side = 4)
mtext("time [week]", side = 4, line = 2)
dev.off()




