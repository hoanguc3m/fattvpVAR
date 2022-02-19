library(fattvpVAR)
numCores = 16
setwd("/Backup/Ongoing/WP11/")
# setwd("/home/hoanguc3m/Downloads/WP11")
##########################################################################

load("G000.RData")
G000_ML_Ent <- MLEnt_TVPGSV(Chain = G000_obj, numCores = numCores)
save(G000_ML_Ent, file = "G000_ML_Ent.RData")
rm(G000_ML_Ent, G000_obj)

load("G001.RData")
G001_ML_Ent <- Chib.LML.TVPGaussian(Chain = G001_obj, numCores = numCores)
save(G000_ML_Ent, file = "G000_ML_Ent.RData")
rm(G000_ML_Ent, G000_obj)


load("G001.RData")
G001_ML_Ent <- MLEnt_TVPGSV(Chain = G001_obj, numCores = numCores)
save(G001_ML_Ent, file = "G001_ML_Ent.RData")
rm(G001_ML_Ent, G001_obj)

load("G010.RData")
G010_ML_Ent <- MLEnt_TVPGSV(Chain = G010_obj, numCores = numCores)
save(G010_ML_Ent, file = "G010_ML_Ent.RData")
rm(G010_ML_Ent, G010_obj)

load("G100.RData")
G100_ML_Ent <- MLEnt_TVPGSV(Chain = G100_obj, numCores = numCores)
save(G100_ML_Ent, file = "G100_ML_Ent.RData")
rm(G100_ML_Ent, G100_obj)

load("G011.RData")
G011_ML_Ent <- MLEnt_TVPGSV(Chain = G011_obj, numCores = numCores)
save(G011_ML_Ent, file = "G011_ML_Ent.RData")

load("G101.RData")
G101_ML_Ent <- MLEnt_TVPGSV(Chain = G101_obj, numCores = numCores)
save(G101_ML_Ent, file = "G101_ML_Ent.RData")

load("G110.RData")
G110_ML_Ent <- MLEnt_TVPGSV(Chain = G110_obj, numCores = numCores)
save(G110_ML_Ent, file = "G110_ML_Ent.RData")

load("G111.RData")
G111_ML_Ent <- MLEnt_TVPGSV(Chain = G111_obj, numCores = numCores)
save(G111_ML_Ent, file = "G111_ML_Ent.RData")

##########################################################################

load("T000.RData")
T000_ML_Ent <- MLEnt_TVPSSV(Chain = T000_obj, numCores = numCores)
save(T000_ML_Ent, file = "T000_ML_Ent.RData")
rm(T000_ML_Ent, T000_obj)


load("T001.RData")
T001_ML_Ent <- MLEnt_TVPSSV(Chain = T001_obj, numCores = numCores)
save(T001_ML_Ent, file = "T001_ML_Ent.RData")
rm(T001_ML_Ent, T001_obj)

load("T010.RData")
T010_ML_Ent <- MLEnt_TVPSSV(Chain = T010_obj, numCores = numCores)
save(T010_ML_Ent, file = "T010_ML_Ent.RData")
rm(T010_ML_Ent, T010_obj)

load("T100.RData")
T100_ML_Ent <- MLEnt_TVPSSV(Chain = T100_obj, numCores = numCores)
save(T100_ML_Ent, file = "T100_ML_Ent.RData")
rm(T100_ML_Ent, T100_obj)

load("T011.RData")
T011_ML_Ent <- MLEnt_TVPSSV(Chain = T011_obj, numCores = numCores)
save(T011_ML_Ent, file = "T011_ML_Ent.RData")

load("T101.RData")
T101_ML_Ent <- MLEnt_TVPSSV(Chain = T101_obj, numCores = numCores)
save(T101_ML_Ent, file = "T101_ML_Ent.RData")

load("T110.RData")
T110_ML_Ent <- MLEnt_TVPSSV(Chain = T110_obj, numCores = numCores)
save(T110_ML_Ent, file = "T110_ML_Ent.RData")

load("T111.RData")
T111_ML_Ent <- MLEnt_TVPSSV(Chain = T111_obj, numCores = numCores)
save(T111_ML_Ent, file = "T111_ML_Ent.RData")


setwd("/home/hoanguc3m/Downloads/WP11/ML")
load("G000_ML_Ent.RData"); load("G001_ML_Ent.RData"); load("G010_ML_Ent.RData"); load("G100_ML_Ent.RData")
load("G011_ML_Ent.RData"); load("G101_ML_Ent.RData"); load("G110_ML_Ent.RData"); load("G111_ML_Ent.RData")

load("T000_ML_Ent.RData"); load("T001_ML_Ent.RData"); load("T010_ML_Ent.RData"); load("T100_ML_Ent.RData")
load("T011_ML_Ent.RData"); load("T101_ML_Ent.RData"); load("T110_ML_Ent.RData"); load("T111_ML_Ent.RData")


xtable::xtable(
  round(rbind(c(G000_ML_Ent$LL, G100_ML_Ent$LL, G010_ML_Ent$LL, G001_ML_Ent$LL, G110_ML_Ent$LL, G101_ML_Ent$LL, G011_ML_Ent$LL, G111_ML_Ent$LL),
              c(T000_ML_Ent$LL, T100_ML_Ent$LL, T010_ML_Ent$LL, T001_ML_Ent$LL, T110_ML_Ent$LL, T101_ML_Ent$LL, T011_ML_Ent$LL, T111_ML_Ent$LL)), digits = 2)
)

