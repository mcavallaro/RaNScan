source('R/Init2.R')
source('R/Evaluate.R')
source("R/plotBaseMap.r")
emmtype = '94.0'

case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"

ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'SAMPLE_DT_numeric')
load(paste0("Data/", emmtype, "_obs.Rdata"))
baseline.matrix = ranScanCreateBaselineMatrix(case.file)
emmtype.factor = ranScanEmmtypeFactor(case.file)
postcode2coord = ranScanPostcodeMap(observation.matrix)

n_cyl = 20000
week.range = c(0, 293)

set.seed(1)
cylinders.94 = ranScanCreateCylinders(observation.matrix,
                                   baseline.matrix * emmtype.factor[[emmtype]], emmtype, week.range,
                                   n_cyl, postcode2coord,  size_factor = 1.1)
set.seed(1)
case.df.tmp.94 = ranScanEvaluate(case.file, cylinders.94, emmtype, p.val.threshold=0.05)

idx = ((case.df.tmp.94$emmtype == emmtype) & !is.na(case.df.tmp.94$x))
case.df.tmp1 = case.df.tmp.94[idx, ]
set.seed(1); tmp.94 = ranScanCluster(case.df.tmp1, emmtype)

png(paste('Fig/', emmtype, '2.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.94, case.df.tmp1, emmtype, legend.pos="bottomright")
dev.off()

png(paste('Fig/', emmtype, 'thresh2.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.94, case.df.tmp1, emmtype, threshold=0.95, legend.pos="bottomright")
dev.off()

png(paste('Fig/', emmtype, 'regions2.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster3(tmp.94, case.df.tmp1, emmtype,
                    legend.pos="topleft",
                    xlim=c(-50, 25), threshold = 0.95)
dev.off()


idxw = case.df.tmp1$warning.score >= 0.95
palette = plasma(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)+1 - min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric), begin = 0, end = 0.75)

png('Fig/_94.0_clust2.png', res=500, width = 5, height = 4, units='in')

par(fig=c(0, 1, 0, 1))
par(mar=c(0, 0, 0, 4))

case.df.tmp1[idxw,]$latitude %>% range
case.df.tmp1[idxw,]$longitude %>% range

plotBaseMap(add=F, xlim=c(-6,0), ylim=c(50.5, 51)) #emm type 94

colors = palette[case.df.tmp1[idxw,]$SAMPLE_DT_numeric - min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric) + 1]
points(latitude ~ longitude,
       data=case.df.tmp1[idxw,],
       col=colors, cex=1.8)

par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
par(mar=c(0,0,0,0))
min_week = min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)

plot(c(0, 1), c(min_week, length(palette)+min_week), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
rasterImage(Cb[seq(length(palette), 1, -1),], 0, min_week, 1, length(palette) + min_week, interpolate = F)
rect(0, min_week, 1, min_week+length(palette))

rr = range(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)
ticks = 119:122
# axis(side = 4, at=ticks, tck=-0.5, labels = NA)
# ticks = ticks[round(seq(1,length(ticks), length.out = 3))]
axis(side = 4, at=ticks, labels=week2Date(ticks))

dev.off()
