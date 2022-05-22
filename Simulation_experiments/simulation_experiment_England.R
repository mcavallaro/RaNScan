## ----setup-----------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Documents/outbreak-detection/")
getwd()


## --------------------------------------------------------------------------------------------------------------------------------------
source('R/utils.R')
source("R/Evaluate.r")


## --------------------------------------------------------------------------------------------------------------------------------------
n.endemic_cases = 5000
size_factor_epi = 4 # multiply by this factor to scale the number of epidemic cases.
label.str = paste(as.character(size_factor_epi), as.character(n.endemic_cases),  sep = '_')
load("Data/population_of_england.RData_")
load('Data/time.factor.RData_')
head(population)

# Set the random generator seed
set.seed(1)
# subsample 10000 postcodes for speed
idx=sort(sample(length(population), 10000))
sample.population = population[idx]
end.matrix = Simulate(sample.population, time.factor[-1][1:100], n.endemic_cases)


## --------------------------------------------------------------------------------------------------------------------------------------
b.matrix = sample.population %o% time.factor[-1][1:100]
b.matrix = b.matrix / sum(b.matrix) * n.endemic_cases


## --------------------------------------------------------------------------------------------------------------------------------------
#geo.location = t(sapply(names(population), postcode.to.location3))
#save(geo.location, file = "Data/geo.location_of_england.RData")
load("Data/geo.location_of_england.RData_")
df.population = cbind(geo.location[idx,], sample.population)
df.cases = cbind(geo.location[idx,], rowSums(end.matrix))
colnames(df.cases) = c('latitude', 'longitude', 'n.cases')
df.cases = as.data.frame(df.cases)
head(df.cases)
df = data.frame(longitude = rep(df.cases[,'longitude'], df.cases[,'n.cases']),
             latitude = rep(df.cases[,'latitude'], df.cases[,'n.cases']))


## --------------------------------------------------------------------------------------------------------------------------------------
idx1 = grepl('AL1', rownames(end.matrix), fixed=T)
idx1 = idx1 | grepl('AL2', rownames(end.matrix), fixed=T)
cat("the number of postcodes starting with AL1 or AL2 is", sum(idx1))
 
# set.seed(1)
save(".Random.seed",file="random_state_seed_1.RData") ## save current RNG state
load("random_state_seed_1.RData")
epi = rpois(n = sum(idx1) * 20, lambda=rep(dnorm(-9:10, sd = 4) * size_factor_epi, rep(sum(idx1), 20)))
save(epi, file = "Data/simulated_larger_epidemic.RData")
cat("the numer of epidemic cases is", sum(epi))

epi.matrix = matrix(epi, ncol = 20)
rownames(epi.matrix) = row.names(end.matrix)[1:sum(idx1)]
colnames(epi.matrix) = 40:59 
epi.matrix[1:3,1:3]
#png(paste0('Fig/epi_time_', label.str, '.png'), width = 4 * 1.2, height = 3, units = 'in', res=400, pointsize = 13)
par(mar=c(4,4,1,1)+0.1)
plot(40:59, colSums(epi.matrix), col=ifelse(colSums(epi.matrix) > 0,'#d62728','white'), xlab='t', pch=20, ylab='no. epidemic cases')
#dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------
all.matrix = end.matrix
all.matrix[1:sum(idx1), as.character(40:59)] = all.matrix[1:sum(idx1), as.character(40:59)] + epi.matrix

# code that requires UK boundaries shape files is commented
# source('R/plotBaseMap.r')
# png(paste0('Fig/end_epi_panel_',label.str,'.png'), width = 4 * 1.2 * 2, height = 3 * 1.2 * 2, units = 'in', res=400, pointsize = 13)
par(mfrow = c(2, 2)) #, mar=c(4,4,1,1)+0.1)
plot(colSums(end.matrix), xlab='t', col='#1f77b4', pch=20, ylab='no. cases', main='endemic')
lines(time.factor[-1][1:100] * n.endemic_cases / sum(time.factor[-1][1:100]), col='#1f77b4')
# mtext('A', side=3, padj=2, at=95, cex=2)
plot(colSums(all.matrix), pch=20, xlab='t', main='endemic + epidemic',  ylab='no. cases')
lines(time.factor[-1][1:100] * n.endemic_cases / sum(time.factor[-1][1:100]), col='#1f77b4')
# mtext('B', side=3, padj=2, at=95, cex=2)

df.cases2 = cbind(geo.location[idx,], rowSums(all.matrix))
df.cases2 = as.data.frame(df.cases2)
df.cases2 = cbind(rownames(geo.location[idx,]), df.cases2)
colnames(df.cases2) = c('Postcode','latitude', 'longitude', 'n.cases')
df2 = data.frame(longitude = rep(df.cases2[,'longitude'], df.cases2[,'n.cases']),
             latitude = rep(df.cases2[,'latitude'], df.cases2[,'n.cases']))

plot(df$longitude, df$latitude,axes=F, pch=20, col='#1f77b4', cex=0.1, xlab=NA, ylab=NA, main='endemic', ylim = c(50,56))
# plotBaseMap(add=T, Wales=F)
# mtext('C', side=3, padj=2, at=1, cex=2)

plot(df2$longitude, df2$latitude,axes=F, pch=20, cex=0.1, xlab=NA, ylab=NA, main='endemic + epidemic', ylim = c(50,56))
# plotBaseMap(add=T, Wales=F)
# mtext('D', side=3, padj=2, at=1, cex=2)


## --------------------------------------------------------------------------------------------------------------------------------------
case.df = as.data.frame(which(all.matrix == 1, arr.ind = TRUE))
case.df$postcode = rownames(all.matrix)[case.df$row]
for (i in 2:max(c(all.matrix))){
  case.df.tmp = as.data.frame(which(all.matrix == i, arr.ind = TRUE))
  case.df.tmp$postcode = rownames(all.matrix)[case.df.tmp$row]
  case.df = rbind(case.df, case.df.tmp[rep(seq_len(NROW(case.df.tmp)), i),])
}

case.df = cbind(case.df, df.population[case.df[,'row'],])
case.df$col = case.df$col - 1
names(case.df) = c('row','SAMPLE_DT_numeric', 'postcode', 'latitude', 'longitude', 'population')
case.df[, c('y','x')] = vlatlong2km(case.df[,c('latitude', 'longitude')])
tail(case.df)


## --------------------------------------------------------------------------------------------------------------------------------------
case.df.epidemic =  which(epi.matrix > 0,  arr.ind = T)
case.df$true_positive = apply(case.df, 1, function(x){ ( as.numeric(x['row']) %in% case.df.epidemic[,'row']) & (as.numeric(x['SAMPLE_DT_numeric']) %in% (case.df.epidemic[,'col'] + 39) )})
# add 39 because epi.matrix starts from time index 40.


## ----echo=TRUE-------------------------------------------------------------------------------------------------------------------------
# set.seed(1)
save(".Random.seed",file="random_state_seed_2.RData") ## save current RNG state
load("random_state_seed_2.RData")
cylinders = CreateCylinders(observation.matrix = all.matrix, baseline.matrix = b.matrix, emmtype = 'sim', week.range = c(0,99), n.cylinders = 1000000, coord.df=df.cases2)
case.df[,c('warning.score', 'low','upp','p.value')] = t(apply(case.df, 1, FUN=warning.score, cylinders))


## --------------------------------------------------------------------------------------------------------------------------------------
head(case.df)


## --------------------------------------------------------------------------------------------------------------------------------------
h = cylinders$t.upp - cylinders$t.low
vol = cylinders$h * cylinders$rho^2
cat("the volume of each cylinder is:", vol[1])
cat("the number of case expected under the baseline model in the cylinders is", quantile(rnorm(100), c(0.05,0.5,0.95)))


## ----echo=TRUE-------------------------------------------------------------------------------------------------------------------------
#set.seed(1)
load("random_state_seed_2.RData")
cylinders2 = CreateCylinders(observation.matrix = all.matrix, baseline.matrix =  b.matrix, emmtype = 'sim', week.range = c(0,99), n.cylinders = 1000000, coord.df=df.cases2,  size_factor = 1.5)
case.df[,c('warning.score2','low2','upp2','p.value2')] = t(apply(case.df, 1, FUN=warning.score, cylinders2))


## ----echo=TRUE-------------------------------------------------------------------------------------------------------------------------
cylinders_a = CreateCylinders(observation.matrix = all.matrix, baseline.matrix =  b.matrix, emmtype = 'sim', week.range = c(0,99), n.cylinders = 500000, coord.df=df.cases2,  size_factor = 1.2)
ws_a = as.data.frame(t(apply(case.df, 1, FUN=warning.score, cylinders_a)))
cylinders_b = CreateCylinders(observation.matrix = all.matrix, baseline.matrix =  b.matrix, emmtype = 'sim', week.range = c(0,99), n.cylinders = 500000, coord.df=df.cases2,  size_factor = 1.5)
ws_b = as.data.frame(t(apply(case.df, 1, FUN=warning.score, cylinders_b)))
names(ws_b) = c('warning.score', 'low','upp','p.value')
names(ws_a) = c('warning.score', 'low','upp','p.value')


## --------------------------------------------------------------------------------------------------------------------------------------
# png(paste0('widths_plot', label.str, '.png'), width = 4.25, height = 3.25, units = 'in', res=600, pointsize = 10)
# par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(density(case.df$upp2 - case.df$low2), main=NA, col='#7f7f7f', xlab='Width of warning-score CI')
lines(density(case.df$upp - case.df$low), col='#1f77b4')
lines(density(ws_b$upp - ws_b$low), col='#7f7f7f', lty=2)
lines(density(ws_a$upp - ws_a$low), col='#1f77b4', lty=2)
#
legend('topright', legend = c("V1, N=10^6", "V2, N=10^6","V1, N=5x10^5", "V2, N=5x10^5"), col=c('#1f77b4', '#7f7f7f', '#1f77b4', '#7f7f7f'), lty=c(1,1,2,2))
# dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------
#case.df.warning.score.old = case.df$warning.score
ctest = cor.test( case.df$warning.score2, case.df$warning.score)
ctest
print(ctest$estimate)
print(ctest$p.value)
colors = ifelse(case.df$true_positive, "#d6272844", "#1f77b444")
# plot(case.df$warning.score, case.df$warning.score2,
#      xlab = 'replicate 1', ylab = 'replicate 2', col=colors, pch=ifelse(case.df$true_positive, 20,1))
#png(paste0('Fig/corrplot', label.str, '.png'), width = 3.25, height = 3.25, units = 'in', res=400, pointsize = 13)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(case.df$warning.score, case.df$warning.score2,
     xlab = expression(paste(italic(w),' (rep. 1)')), ylab = expression(paste(italic(w)," (rep. 2)")), col=colors, pch=ifelse(case.df$true_positive, 20,1))
legend('bottomright', c('epidemic','endemic'), col=c("#d6272855", "#1f77b455"), pch=c(20,1))
#dev.off()


## ----message=TRUE, warning=FALSE-------------------------------------------------------------------------------------------------------
library(pROC)
ROC = roc(case.df$true_positive, case.df$warning.score)
print(ROC)
plot(ROC)
L = NROW(case.df)
col = '#00000033'
#png(paste0('Fig/ROC_',label.str,'.png'), width = 3.25, height = 3.25, units = 'in', res=400, pointsize = 13)
par(mfrow=c(1,1), mar=c(4,4,0,0))
auc = numeric(0)
#plot(ROC, col='black', identity.col='black',  grid.col='red')
for (i in 1:5){
  cc = roc(true_positive ~ warning.score, case.df[sample(1:L, L, replace = T),])
  if(i < 11){
#    plot(cc, add=T, col=col, identity.col='black')    
  }
  auc = c(auc, as.numeric(cc['auc'][[1]]))
}
#dev.off()
print(quantile(auc, c(0.025,0.5,0.975) ))


## --------------------------------------------------------------------------------------------------------------------------------------
# png(paste0('Fig/ws_time_',label.str,'.png'), width = 4 * 1.2, height = 3 * 1.2, units = 'in', res=400, pointsize = 13)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(case.df$SAMPLE_DT_numeric, case.df$warning.score, xlab='t',
     ylab=expression(italic(w)),
     col=ifelse(case.df$true_positive, "#d6272855", "#1f77b455"),
     pch=ifelse(case.df$true_positive, 20, 1))
legend('bottomright', c('epidemic','endemic'), col=c("#d6272855", "#1f77b455"),
       pch=c(20,1))
# dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------
# palette = colorRampPalette(c('#1f77b4', '#d62728'))(100)
idx1 = (case.df$SAMPLE_DT_numeric > 40 ) & ( case.df$SAMPLE_DT_numeric < 60) & (case.df$warning.score < 0.95)
#colors = palette[as.integer(case.df[idx,]$warning.score * 100)+1]
idx2 = (case.df$SAMPLE_DT_numeric > 40 ) & ( case.df$SAMPLE_DT_numeric < 60) & (case.df$warning.score > 0.95)
# plot(case.df[idx1,]$longitude, case.df[idx1,]$latitude, col='#1f77b4',axes=F, pch=1, cex=1, ylab = NA, xlab=NA)
# points(case.df[idx2,]$longitude, case.df[idx2,]$latitude, col='#d62728', pch=20, cex=1, ylab = NA, xlab=NA)
#
# png(paste0('Fig/satscan_',label.str,'.png'), width =  4 *1.2 / 1.5, height = 3 *1.2 /1.5, units = 'in', res=400, pointsize = 20)
par(mfrow = c(1, 1), mar=c(0,0,0,0) + 0.5)
idx3=case.df$warning.score > 0.95
plot(case.df[!idx3,]$longitude, case.df[!idx3,]$latitude, col='#1f77b4',axes=F, pch=1, cex=1, ylab = NA, xlab=NA,  xlim=c(-0.58, 0.01), ylim=c(51.57, 51.8) )
points(case.df[case.df$true_positive,]$longitude, case.df[case.df$true_positive,]$latitude, col='#d62728', pch=20, cex=1, ylab = NA, xlab=NA)
points(case.df[idx3,]$longitude, case.df[idx3,]$latitude, col='#ff7f0e', pch=4, cex=1, ylab = NA, xlab=NA)

# plotBaseMap(add=T, Wales = F, onlyregion = T)
# sim.col=st_read('shape_files/satscan_output.shp')
# plot(sim.col[1], add=T, col=NA, border='#7f7f7f')
# dev.off()


## --------------------------------------------------------------------------------------------------------------------------------------
# png(paste0('Fig/end_epi_spatial_',label.str,'.png'), width =  4 *1.6, height = 3 *1.6, units = 'in', res=600, pointsize = 13)
par(mfrow = c(1, 1), mar=c(1,1,1,1)+0.1)

plot(case.df[idx1,]$longitude, case.df[idx1,]$latitude, col='#1f77b4',axes=F, pch=1, cex=1, ylab = NA, xlab=NA, ylim = c(50,56), xlim = c(-5.9, 1.8))
points(case.df[idx2,]$longitude, case.df[idx2,]$latitude, col='#ff7f0e', pch=4, cex=1, ylab = NA, xlab=NA)
# plotBaseMap(add=T, Wales = F)
legend('topleft',
       c( expression(italic(w)>=0.95),
          expression(italic(w)<0.95),
          "True outbreak"), col=c("#ff7f0e", "#1f77b4", '#d62728'), pch=c(4,1,20), box.col = "white")
# dev.off()


## ----echo=TRUE-------------------------------------------------------------------------------------------------------------------------
save(".Random.seed",file="random_state_seed_3.RData")
load("random_state_seed_3.RData")
for (week in 40:59){
  week.range = as.character(0:week)
  cylinders = CreateCylinders(observation.matrix = all.matrix[,week.range], baseline.matrix = b.matrix[,week.range], emmtype = 'sim', week.range = week.range, n.cylinders = 1000000, coord.df=df.cases2,  size_factor = 1.2)
  case.df[,paste0('warning.score', as.character(week))] = apply(case.df, 1, FUN=warning.score, cylinders)[1,]
}


## --------------------------------------------------------------------------------------------------------------------------------------
library(viridisLite)
palette=inferno(length(40:59),begin=0.15, end=0.85, direction=-1)

plot(c(40,59), range(case.df[case.df$true_positive,11:30]), col='white', xlab = expression(tau),
     ylab=expression(italic(w)), axes = FALSE)
#abline(h=0.95, col='#7f7f7f')

for(i in 1:NROW(case.df[case.df$true_positive,  ])){
  t1 = case.df[case.df$true_positive,][i,'SAMPLE_DT_numeric']
  id = t1 - 40 + 18
  # 40 is the week of the start of the outbreak
  # 18 is the number of columns before the column `warning.score40`
  color = palette[t1-39]
  if (t1 <= 59){
  #  print(c(dim(case.df[case.df$true_positive,]),  i, id:31))
#    print( unlist(case.df[case.df$true_positive,][i, id:31]))
    lines(t1:59, case.df[case.df$true_positive,][i, id:37], col=color)
    points(t1:59, case.df[case.df$true_positive,][i, id:37], pch=20, cex=0.4, col=color)
  }
}
axis(2, at=c(0,0.25,0.5,0.75,1))
axis(1, at = seq(40,59))

