
emmtype = '108.1'

# source("ranScanInit.r")
case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"
load(paste0("~/Outbreak/Data/", emmtype, "_obs.Rdata"))
baseline.matrix = ranScanCreateBaselineMatrix(case.file)
emmtype.factor = ranScanEmmtypeFactor(case.file)

n_cyl = 20000
week.range = c(1, 298)

set.seed(1)
cylinders.1081 = ranScanCreateCylinders(observation.matrix,
                                   baseline.matrix * emmtype.factor[[emmtype]], emmtype, week.range,
                                   n_cyl, postcode2coord,  size_factor = sizefactor)
set.seed(1)
case.df.tmp.1081 = ranScanEvaluate(case.file, cylinders.1081, emmtype, p.val.threshold=0.05)

idx = ((case.df.tmp.1081$emmtype == emmtype) & !is.na(case.df.tmp.1081$x))
case.df.tmp1 = case.df.tmp.1081[idx, ]

set.seed(1); tmp.1081 = ranScanCluster(case.df.tmp1, emmtype)

png(paste('Manuscript/', emmtype, 'regions.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster3(tmp.1081, case.df.tmp1, emmtype, Legend = T, xlim=c(-36,15), ylim=c(-15,17))
dev.off()

png(paste('Manuscript/', emmtype, '.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.1081, case.df.tmp1, emmtype)
dev.off()

png(paste('Manuscript/', emmtype, 'thresh.png', sep='_'),
    res=600, width=5, height=4, units="in")
ranScanPlotCluster(tmp.1081, case.df.tmp1, emmtype, threshold=0.95)
dev.off()

# idxx = (tmp[,1] < 9.5) & (tmp[,1] > 3.5)
# idxy = (tmp[,2] < -3.5) & (tmp[,2] > -9.5)

# idxx1 = (tmp[,1] < 5) & (tmp[,1] > 3)
# idxy1 = (tmp[,2] < 3.5) & (tmp[,2] > 0)
# 
# idxx2 = (tmp[,1] < 10) & (tmp[,1] >4)
# idxy2 = (tmp[,2] < 10) & (tmp[,2] > 3.5)
# idxx = idxx1 | idxx2
# idxy = idxy1 | idxy2
idxw = case.df.tmp1$warning.score > 0.95
# idxw1 = case.df.tmp1[idxx1 & idxy1,]$warning.score > 0.95
# idxw2 = case.df.tmp1[idxx2 & idxy2,]$warning.score > 0.95
palette = plasma(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric) -
                    min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric), begin = 0, end = 0.75)

print(c(max(case.df.tmp1[idxw,]$SAMPLE_DT_numeric),
        min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)))

png('Manuscript/_108.1clust.png', res=500, width = 5, height = 4, units='in')

par(fig=c(0, 1, 0, 1))
par(mar=c(0, 0, 0, 4))


y.range = range(case.df.tmp1[idxw,]$latitude)
x.range = range(case.df.tmp1[idxw,]$longitude)

plotBaseMap(add=F, xlim=x.range, ylim=y.range) 

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

min_week = min(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)

plot(c(0, 1), c(min_week, len(palette)+min_week), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
rasterImage(Cb[seq(length(palette), 1, -1),], 0, min_week, 1, len(palette) + min_week)
rect(0, min_week, 1, min_week+length(palette))

rr = range(case.df.tmp1[idxw,]$SAMPLE_DT_numeric)
ticks = seq(rr[1], rr[2])

axis(side = 4, at=ticks, labels=week2Date(ticks))


dev.off()



