library(fattvpVAR)
numCores = 16
setwd("/home/hoanguc3m/Downloads/WP11/")

##########################################################################

load("G000.RData")
G000_ML <- ML_GaussTVPSV(Chain = G000_obj, numCores = numCores)
save(G000_ML, file = "G000_ML.RData")
rm(G000_ML, G000_obj)


load("G001.RData")
G001_ML <- ML_GaussTVPSV(Chain = G001_obj, numCores = numCores)
save(G001_ML, file = "G001_ML.RData")
rm(G001_ML, G001_obj)

load("G010.RData")
G010_ML <- ML_GaussTVPSV(Chain = G010_obj, numCores = numCores)
save(G010_ML, file = "G010_ML.RData")
rm(G010_ML, G010_obj)

load("G100.RData")
G100_ML <- ML_GaussTVPSV(Chain = G100_obj, numCores = numCores)
save(G100_ML, file = "G100_ML.RData")
rm(G100_ML, G100_obj)

load("G011.RData")
G011_ML <- ML_GaussTVPSV(Chain = G011_obj, numCores = numCores)
save(G011_ML, file = "G011_ML.RData")

load("G101.RData")
G101_ML <- ML_GaussTVPSV(Chain = G101_obj, numCores = numCores)
save(G101_ML, file = "G101_ML.RData")

load("G110.RData")
G110_ML <- ML_GaussTVPSV(Chain = G110_obj, numCores = numCores)
save(G110_ML, file = "G110_ML.RData")

load("G111.RData")
G111_ML <- ML_GaussTVPSV(Chain = G111_obj, numCores = numCores)
save(G111_ML, file = "G111_ML.RData")

##########################################################################

load("T000.RData")
T000_ML <- ML_StudentTVPSV(Chain = T000_obj, numCores = numCores)
save(T000_ML, file = "T000_ML.RData")
rm(T000_ML, T000_obj)


load("T001.RData")
T001_ML <- ML_StudentTVPSV(Chain = T001_obj, numCores = numCores)
save(T001_ML, file = "T001_ML.RData")
rm(T001_ML, T001_obj)

load("T010.RData")
T010_ML <- ML_StudentTVPSV(Chain = T010_obj, numCores = numCores)
save(T010_ML, file = "T010_ML.RData")
rm(T010_ML, T010_obj)

load("T100.RData")
T100_ML <- ML_StudentTVPSV(Chain = T100_obj, numCores = numCores)
save(T100_ML, file = "T100_ML.RData")
rm(T100_ML, T100_obj)

load("T011.RData")
T011_ML <- ML_StudentTVPSV(Chain = T011_obj, numCores = numCores)
save(T011_ML, file = "T011_ML.RData")

load("T101.RData")
T101_ML <- ML_StudentTVPSV(Chain = T101_obj, numCores = numCores)
save(T101_ML, file = "T101_ML.RData")

load("T110.RData")
T110_ML <- ML_StudentTVPSV(Chain = T110_obj, numCores = numCores)
save(T110_ML, file = "T110_ML.RData")

load("T111.RData")
T111_ML <- ML_StudentTVPSV(Chain = T111_obj, numCores = numCores)
save(T111_ML, file = "T111_ML.RData")


xtable::xtable(
  round(rbind(c(G000_ML$LL, G100_ML$LL, G010_ML$LL, G001_ML$LL, G110_ML$LL, G101_ML$LL, G011_ML$LL, G111_ML$LL),
              c(T000_ML$LL, T100_ML$LL, T010_ML$LL, T001_ML$LL, T110_ML$LL, T101_ML$LL, T011_ML$LL, T111_ML$LL)), digits = 2)
)

