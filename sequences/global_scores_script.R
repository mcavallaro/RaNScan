
# INITIALISE:
source("ranScanInit2.r")
case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"
# 
case.df = ranScanInit(case.file)
emmtypes = c("1.6", "11.1", "12.87", "25.1", "5.132", "169.3", "5", "6.59","102.2",
                "3.119","49.0","1.29",
                "18.0","8.0","5.5", "1.25","90.2", "g245.0", "12.37", "58.0", "5.129","g652.0","76.0","73.0","1.3","121.0","80.0",
                "83.13", "92.0", "9.0", "6.4", "33.0","22.0", "g485.0","81.0", "5.3", "2.0",
                "5.23", "82.0", "44.0", "11.0", "77.0", "6.0", "3.1", "75.0",
                "108.1", "87.0", "3.93", "66.0", "4.0", "94.0", "28.0", "12.0",
                "89.0", "1.0")
ranScanCreateObservationMatrices(case.file, emmtypes)
load("Data/1.0_obs.Rdata", verbose=1)
postcode2coord = ranScanPostcodeMap(observation.matrix)
baseline.matrix = ranScanCreateBaselineMatrix(case.file)
emmtype.factor.all = ranScanEmmtypeFactor(case.file)

time.factor = ranScanTimeFactor(case.file)
# #
# sum(baseline.matrix * emmtype.factor[['1.0']], na.rm=T)
# sum(observation.matrix, na.rm=T)
# 
# load("Data/12.0_obs.Rdata")
# sum(baseline.matrix * emmtype.factor[['12.0']], na.rm=T)
# sum(observation.matrix, na.rm=T)
# source("ranScanEvaluate.r")
# ranScanCreateObservationMatrices(case.file, c('33.0', '108.1', '12.0'))
# ranScanCreateObservationMatrices(case.file, c('44.0'))

source("ranScanEvaluate.r")


global.scores = vector(mode = 'numeric', length = length(emmtypes))
sizefactor = 1.2
week.range = c(150, 293)
for (i in 1:length(emmtypes)){
  emmtype = emmtypes[i]
  t1 = Sys.time()
  load(paste0("~/Outbreak/Data/", emmtype, "_obs.Rdata"))
  if (!(emmtype == attributes(observation.matrix)$emmtype )){
    stop(emmtype)
  }
  n_cyl = min(10000, sum(observation.matrix) * 10)
  cylinders = ranScanCreateCylinders(observation.matrix,
                                     baseline.matrix * emmtype.factor.all[[emmtype]],
                                     emmtype, week.range,
                                     n_cyl, postcode2coord, size_factor = sizefactor)
  
  case.df.tmp = ranScanEvaluate(case.file, cylinders, emmtype, p.val.threshold=0.05)
  global.score = mean(case.df.tmp[case.df.tmp$emmtype == emmtype, ]$warning.score)
  writeLines(sprintf("The global score for emm type %s is %f", emmtype, global.score))
  global.scores[i] = global.score
  print(Sys.time() - t1)
}


global.scores = data.frame(emmtype = as.character(emmtypes), global.score = global.scores)


global.scores$emmtype = as.character(global.scores$emmtype)
global.scores$prevalence = apply(global.scores, 1, function(x){sum(case.df$emmtype == x['emmtype'])}) / sum(postcode2coord$Total, na.rm=T)
global.scores$fraction = global.scores$prevalence / sum(global.scores$prevalence)

idx=order(global.scores$global.score, decreasing = T)

write.table(global.scores[idx,], file = "~/Outbreak/sequences/global_scores.csv", quote = F, sep = ",", row.names = F, col.names = TRUE)



