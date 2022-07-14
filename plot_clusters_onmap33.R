source('R/Init2.R')
source('R/Evaluate.R')
source("R/plotBaseMap.r")
emmtype = '33.0'

case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"
ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'SAMPLE_DT_numeric')
#ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'RECEPT_DT_numeric')
load(paste0("Data/", emmtype, "_obs.Rdata"))
baseline.matrix = ranScanCreateBaselineMatrix(case.file)
emmtype.factor = ranScanEmmtypeFactor(case.file)

postcode2coord = ranScanPostcodeMap(observation.matrix)

sizefactor=1.2
n_cyl = 20000
week.range = c(0, 293)

set.seed(1)
cylinders.33 = ranScanCreateCylinders(observation.matrix,
                                   baseline.matrix * emmtype.factor[[emmtype]], emmtype, week.range,
                                   n_cyl, postcode2coord,  size_factor = sizefactor)
set.seed(1)
case.df.tmp.33 = ranScanEvaluate(case.file, cylinders.33, emmtype, p.val.threshold=0.05)

idx = ((case.df.tmp.33$emmtype == emmtype) & !is.na(case.df.tmp.33$x))
case.df.tmp1 = case.df.tmp.33[idx, ]

set.seed(1); tmp.33 = ranScanCluster(case.df.tmp1, emmtype)

png(paste('Fig/', emmtype, 'regions2.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster3(tmp.33, case.df.tmp1, emmtype, xlim=c(-2200, 1500))
dev.off()

png(paste('Fig/', emmtype, '2.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.33, case.df.tmp1, emmtype)
dev.off()

png(paste('Fig/', emmtype, 'thresh2.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.33, case.df.tmp1, emmtype, threshold=0.95)
dev.off()

# # idxx = (tmp[,1] < 9.5) & (tmp[,1] > 3.5)
# # idxy = (tmp[,2] < -3.5) & (tmp[,2] > -9.5)
# 
# idxx1 = (tmp[,1] < 5) & (tmp[,1] > 3)
# idxy1 = (tmp[,2] < 3.5) & (tmp[,2] > 0)
# idxx2 = (tmp[,1] < 10) & (tmp[,1] >4)
# idxy2 = (tmp[,2] < 10) & (tmp[,2] > 3.5)
# idxx = idxx1 | idxx2
# idxy = idxy1 | idxy2

idxw = case.df.tmp1$warning.score > 0.95

palette = plasma(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric) -
                   min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric), begin = 0, end = 0.75)

print(c(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric), min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)))

case.df.tmp1[idxw,]$latitude %>% range
case.df.tmp1[idxw,]$longitude %>% range
# plotBaseMap(add=F, xlim=c(-5.5,1.5), ylim=c(50,52)) #emm type 94

png('Fig/_33.0_clust2.png', res=500, width = 5, height = 4, units='in')
par(fig=c(0, 1, 0, 1))
par(mar=c(0, 0, 0, 4))

min_week = min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)

plotBaseMap(add=F, xlim=c(-2.2, -0.3), ylim=c(51.8,54), onlyregion = F)

colors = palette[case.df.tmp1[idxw,]$SAMPLE_DT_numeric - min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)]
points(latitude ~ longitude,
       data=case.df.tmp1[idxw,],
       col=colors, cex=1.7)

par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
par(mar=c(0,0,0,0))
plot(c(0, 1), c(min_week, length(palette) + min_week), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')

Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
rasterImage(Cb[seq(length(palette), 1, -1),], 0, min_week, 1, length(palette) + min_week)
rect(0, min_week, 1, min_week+length(palette))

rr = range(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)
ticks = seq(rr[1], rr[2])
axis(side = 4, at=ticks, tck=-0.5, labels = NA)
ticks = ticks[round(seq(1,length(ticks), length.out = 3))]
axis(side = 4, at=ticks, tck=-1.1, labels=week2Date(ticks))

dev.off()
