## ----warning=TRUE, include=FALSE-----------------------------------------------------------------------------
source("ranScanInit2.r")
case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx"
case.df = ranScanInit(case.file)
print(head( unique(case.df$`Patient Postcode`), 13 ))

## ------------------------------------------------------------------------------------------------------------
emmtype='44.0'

ranScanCreateObservationMatrices(case.file, emmtype, date.time.field = 'RECEPT_DT_numeric')

load(paste0("~/Outbreak/Data/", emmtype, "_obs.Rdata"))
#time.factor = ranScanTimeFactor(case.file)

print(head(row.names(observation.matrix), 13))

## ------------------------------------------------------------------------------------------------------------
baseline.matrix = ranScanCreateBaselineMatrix(case.file, date.time.field = 'RECEPT_DT_numeric', save.on.disk = FALSE)
print(head(row.names(baseline.matrix),13))


## ------------------------------------------------------------------------------------------------------------
init=Sys.time()
postcode2coord = ranScanPostcodeMap(observation.matrix)
print(Sys.time()-init)
print(as.character(postcode2coord[1:13,1]))

## ------------------------------------------------------------------------------------------------------------
emmtype.factor = ranScanEmmtypeFactor.tau(case.file, emmtype, date.time.field = 'RECEPT_DT_numeric')


## ----include=FALSE-------------------------------------------------------------------------------------------
source("ranScanEvaluate.r")


## ----include=FALSE-------------------------------------------------------------------------------------------
sizefactor = 1.2
n_cyl = 50000 # min(50000, sum(obs.matrix.typed[-13,]) * 10 )

# moltiplicazione columnwise
baseline.matrix.tmp = baseline.matrix * rep(c(0,0,emmtype.factor[[emmtype]]), rep(NROW(baseline.matrix), NCOL(baseline.matrix) ))
postcode2coord = ranScanPostcodeMap(observation.matrix)

weeks = 96:max(as.integer(colnames(observation.matrix)), na.rm = T)
#Warning message: NAs introduced by coercion 

n.rows = sum(case.df$emmtype == emmtype)
n.cols = length(weeks)
to.export = as.data.frame(matrix(NA, nrow = n.rows, ncol = n.cols))
names(to.export) = mapply(function(x){week2Date(x)}, as.character(weeks))
#
for (week in weeks){
  week.min = max(0, week - 150)
  week.range = as.character(week.min:week)
  cat('weeks:', range(as.integer(week.range)),'\t')

  cylinders = ranScanCreateCylinders(observation.matrix[,week.range],
                                     baseline.matrix.tmp[,week.range],
                                     emmtype, week.range, n_cyl, postcode2coord, size_factor = sizefactor)
  
  case.df.tmp = ranScanEvaluate(case.file, cylinders, emmtype, p.val.threshold = 0.05, date.time.field = 'RECEPT_DT_numeric')
  ws = case.df.tmp[case.df.tmp$emmtype == emmtype,]$warning.score
  i = week - weeks[1] + 1
  to.export[, i] = ws
}

## ------------------------------------------------------------------------------------------------------------
dataframe.to.export = cbind(case.df.tmp[case.df.tmp$emmtype == emmtype,], to.export)
dataframe.to.export$warning.score = NULL
write.csv(dataframe.to.export, file = paste0('emmtype_', emmtype, '_I_no_delay_receipt.csv'))

