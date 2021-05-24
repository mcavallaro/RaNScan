#load("~/Outbreak/iGAS.Rdata")

source("utils.R")

dat['emmtypes2']  = apply(dat['emm gene type'],1, FUN=split)

idx = which(dat['emmtypes2']  == "c74a.0\r")
dat$emmtypes2[idx[1]]  = 'c74a.0'
dat$emmtypes2[idx[2]]  = 'c74a.0'

length(unique(dat$emmtypes2) )
#there are 183 unique emmtypes

idx = dat$`Sterile Site Y N` == "#s"
idx2 = !is.na(dat$`Sterile Site Y N`)
dat = dat[idx & idx2, ]

# Insert coordinates
UKpostcodes = read.csv("ukpostcodes.csv", stringsAsFactors = F)


dat[,c("latitude", "longitude")] = t(apply(dat, 1, postcode.to.location, UKpostcodes))


dat$RECEPT_DT_numeric = as.integer(difftime(dat$RECEPT_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
dat$SAMPLE_DT_numeric = as.integer(difftime(dat$SAMPLE_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))

idx = (dat$SAMPLE_DT_numeric > -1) & !is.na(dat$SAMPLE_DT_numeric)
dat<-dat[idx,]

idx = is.na(dat$emmtypes2)
dat$emmtypes2[idx] = 'NA'

idx = is.na(dat$`Patient Postcode`)
dat$`Patient Postcode`[idx] = 'NA'

# idx = is.na(dat$latitude)
# dat$latitude[idx] = 'NA'
# 
# idx = is.na(dat$longitude)
# dat$longitude[idx] = 'NA'


dat$diff = dat$SAMPLE_DT_numeric - dat$RECEPT_DT_numeric

save("dat", file="iGAS_coords_full.Rdata")
