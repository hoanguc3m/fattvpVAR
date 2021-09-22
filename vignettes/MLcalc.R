library(fattvpVAR)
numCores = 4
setwd("/home/hoanguc3m/Downloads/WP11/")
#
#
# load("G000.RData")
# G000_ML <- ML_GaussTVPSV(Chain = G000_obj, numCores = numCores)
# save(G000_ML, file = "G000_ML.RData")
#
#
#
load("G001.RData")
G001_ML <- ML_GaussTVPSV(Chain = G001_obj, numCores = numCores)
save(G001_ML, file = "G001_ML.RData")

load("G010.RData")
G010_ML <- ML_GaussTVPSV(Chain = G010_obj, numCores = numCores)
save(G010_ML, file = "G010_ML.RData")

load("G100.RData")
G100_ML <- ML_GaussTVPSV(Chain = G100_obj, numCores = numCores)
save(G100_ML, file = "G100_ML.RData")

load("G011.RData")
G011_ML <- ML_GaussTVPSV(Chain = G011_obj, numCores = numCores)
save(G011_ML, file = "G011_ML.RData")


# load("G111.RData")
# G111_ML <- ML_GaussTVPSV(Chain = G111_obj, numCores = numCores)
# save(G111_ML, file = "G111_ML.RData")
