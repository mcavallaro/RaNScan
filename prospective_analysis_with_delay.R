## ----include=FALSE-----------------------------------------------------------------------------------------------
source("ranScanInit2.r")
case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"
case.df = ranScanInit(case.file)


## ----------------------------------------------------------------------------------------------------------------
emmtype='44.0' ##44.0 #1.0 #33.0
n.weeks = 150
starting.week = 80
#case.df$`Patient Postcode` %>% unique %>% head(14)


## ----------------------------------------------------------------------------------------------------------------
ranScanCreateObservationMatrices.delay(case.file, emmtype, starting.week = starting.week, n.weeks = n.weeks)


## ----------------------------------------------------------------------------------------------------------------
load(sprintf("Data/%s_typed_obs.Rdata", emmtype), verbose = 1)
load('Data/untyped_obs.Rdata', verbose = 1)

print(head(names(observation.matrices.untyped)))
observation.matrices.untyped[['150']][1:5,1:5]
observation.matrices.untyped[[10]] %>% rownames %>% head
observation.matrices.typed[[emmtype]][['150']][1:5,1:5]
case.df$`Patient Postcode` %>% unique %>% head(14)


## ----------------------------------------------------------------------------------------------------------------
init=Sys.time()
postcode2coord = ranScanPostcodeMap(observation.matrices.untyped[[1]])
print(Sys.time()-init)


## ----------------------------------------------------------------------------------------------------------------
#time.factor = ranScanTimeFactor(case.file)


## ----------------------------------------------------------------------------------------------------------------
baseline.matrix = ranScanCreateBaselineMatrix(case.file, save.on.disk = FALSE)


## ----------------------------------------------------------------------------------------------------------------
emmtype.factor = ranScanEmmtypeFactor.tau(case.file, emmtype)
delay.factor = ranScanEmmtypeFactor.delay_(case.file, starting.week = starting.week, n.weeks = n.weeks)

## ----include=FALSE-----------------------------------------------------------------------------------------------
source("ranScanEvaluate.r")


## ----include=FALSE-----------------------------------------------------------------------------------------------
sizefactor = 1.2
n_cyl = 50000 # min(50000, sum(obs.matrix.typed[-13,]) * 10 )
weeks = attributes(observation.matrices.untyped)$names
n.rows = sum(case.df$emmtype == emmtype)
n.cols = length(weeks)

to.export.strategy1 = as.data.frame(matrix(NA, nrow = n.rows, ncol = n.cols))
to.export.strategy2 = as.data.frame(matrix(NA, nrow = n.rows, ncol = n.cols))
names(to.export.strategy1) = mapply(function(x){week2Date(x)}, as.character(weeks))
names(to.export.strategy2) = mapply(function(x){week2Date(x)}, as.character(weeks))

dly.factor = delay.factor
  
for (week in weeks){
  writeLines(week)
  
  cols = colnames(observation.matrices.untyped[[week]])
  emt.factor = emmtype.factor[[emmtype]][cols]
  
  obs.matrix.typed = observation.matrices.typed[[emmtype]][[week]]
  obs.matrix.untyped = observation.matrices.untyped[[week]] #* rep(emt.factor,
                                                                  # rep(NROW(observation.matrices.untyped[[week]]),
                                                                  #     NCOL(observation.matrices.untyped[[week]])))
  obs.matrix = obs.matrix.typed + obs.matrix.untyped
  
  bsl.matrix.typed = baseline.matrix[,cols] * rep((1 - dly.factor) * emt.factor,
                                                  rep(NROW(baseline.matrix[,cols]),
                                                      NCOL(baseline.matrix[,cols])))
  bsl.matrix.untyped = baseline.matrix[,cols] * rep(dly.factor,
                                                    rep(NROW(baseline.matrix[,cols]),
                                                        NCOL(baseline.matrix[,cols])))
  bsl.matrix = bsl.matrix.typed + bsl.matrix.untyped
  
  week.range = range(as.integer(cols))

  cylinders.strategy1 = ranScanCreateCylinders.delay(obs.matrix.typed, bsl.matrix.typed,
                                           obs.matrix.untyped, bsl.matrix.untyped,
                                           emmtype, week.range,
                                           n_cyl, postcode2coord, size_factor = sizefactor)
  #
  cylinders.strategy2 = ranScanCreateCylinders(obs.matrix, bsl.matrix,
                                     emmtype, week.range,
                                     n_cyl, postcode2coord, size_factor = sizefactor)
  #
  case.df.tmp = ranScanEvaluate(case.file, cylinders.strategy1, emmtype, p.val.threshold=0.05)
  ws = case.df.tmp[case.df.tmp$emmtype == emmtype,]$warning.score
  i = as.integer(week) - as.integer(weeks[1]) + 1
  to.export.strategy1[, i] = ws
  
  case.df.tmp = ranScanEvaluate(case.file, cylinders.strategy2, emmtype, p.val.threshold=0.05)
  ws = case.df.tmp[case.df.tmp$emmtype == emmtype,]$warning.score
  i = as.integer(week) - as.integer(weeks[1]) + 1
  to.export.strategy2[, i] = ws  
  
}


dataframe.to.export = cbind(case.df.tmp[case.df.tmp$emmtype == emmtype,], to.export.strategy1)
dataframe.to.export$warning.score = NULL
write.csv(dataframe.to.export, file = paste0('emmtype_', emmtype, '_2_strategy1_delay.csv'))

dataframe.to.export = cbind(case.df.tmp[case.df.tmp$emmtype == emmtype,], to.export.strategy2)
dataframe.to.export$warning.score = NULL
write.csv(dataframe.to.export, file = paste0('emmtype_', emmtype, '_2_strategy2_delay.csv'))
rm(dataframe.to.export)



