rm(list=ls())
source('utils.R')
postcode.df = read.table("Data/ukpostcodes.csv", stringsAsFactors = F, sep=',', header = T)
Population = read.table("Data/Postcode_Estimates_Table_1.csv", header=T, sep = ',', stringsAsFactors = F)
names(Population)[1] = 'Patient Postcode'
Population[,c("latitude", "longitude")] = t(apply(Population, 1, postcode.to.location, postcode.df))
Population$Patient.Postcode = apply(Population, 1, function(x){gsub(' ', '', x[1], fixed=T)})

write.table(Population, file='Data/Postcode_Estimates_Table_with_coordinates3.csv', quote = F, row.names = F, sep=',')
