metRates  <-  read.csv('data/data.csv', header=TRUE, stringsAsFactors=FALSE)
metRates  <-  metRates[complete.cases(metRates), ]

# correct log values from log 10 to ln
metRates$lnMass  <-  log(10^metRates$logM)
metRates$lnRate  <-  log(10^metRates$logVo2)

# create inverse temperature column
metRates$TempK  <-  metRates$Temp + 273.15
metRates$invKT  <-  1 / 8.62e-5 * (1 / 293.15 - 1 / metRates$TempK)

# factor temperature column for ancova
metRates$Temp         <-  as.factor(metRates$Temp)
