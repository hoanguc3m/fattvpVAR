library(fattvpVAR)
numCores = 4
setwd("/Backup/Ongoing/WP11/")
 setwd("/home/hoanguc3m/Downloads/WP11")

##########################################################################

load("G000.RData")
G000_ML_Chib <- MLChib_TVPGSV(Chain = G000_obj, numCores = numCores)
save(G000_ML_Chib, file = "G000_ML_Chib.RData")
rm(G000_ML_Chib, G000_obj)


load("G001.RData")
G001_ML_Chib <- MLChib_TVPGSV(Chain = G001_obj, numCores = numCores)
save(G001_ML_Chib, file = "G001_ML_Chib.RData")
rm(G001_ML_Chib, G001_obj)

load("G010.RData")
G010_ML_Chib <- MLChib_TVPGSV(Chain = G010_obj, numCores = numCores)
save(G010_ML_Chib, file = "G010_ML_Chib.RData")
rm(G010_ML_Chib, G010_obj)

load("G100.RData")
G100_ML_Chib <- MLChib_TVPGSV(Chain = G100_obj, numCores = numCores)
save(G100_ML_Chib, file = "G100_ML_Chib.RData")
rm(G100_ML_Chib, G100_obj)

load("G011.RData")
G011_ML_Chib <- MLChib_TVPGSV(Chain = G011_obj, numCores = numCores)
save(G011_ML_Chib, file = "G011_ML_Chib.RData")

load("G101.RData")
G101_ML_Chib <- MLChib_TVPGSV(Chain = G101_obj, numCores = numCores)
save(G101_ML_Chib, file = "G101_ML_Chib.RData")

load("G110.RData")
G110_ML_Chib <- MLChib_TVPGSV(Chain = G110_obj, numCores = numCores)
save(G110_ML_Chib, file = "G110_ML_Chib.RData")

load("G111.RData")
G111_ML_Chib <- MLChib_TVPGSV(Chain = G111_obj, numCores = numCores)
save(G111_ML_Chib, file = "G111_ML_Chib.RData")

##########################################################################

load("T000.RData")
T000_ML_Chib <- MLChib_TVPSSV(Chain = T000_obj, numCores = numCores)
save(T000_ML_Chib, file = "T000_ML_Chib.RData")
rm(T000_ML_Chib, T000_obj)


load("T001.RData")
T001_ML_Chib <- MLChib_TVPSSV(Chain = T001_obj, numCores = numCores)
save(T001_ML_Chib, file = "T001_ML_Chib.RData")
rm(T001_ML_Chib, T001_obj)

load("T010.RData")
T010_ML_Chib <- MLChib_TVPSSV(Chain = T010_obj, numCores = numCores)
save(T010_ML_Chib, file = "T010_ML_Chib.RData")
rm(T010_ML_Chib, T010_obj)

load("T100.RData")
T100_ML_Chib <- MLChib_TVPSSV(Chain = T100_obj, numCores = numCores)
save(T100_ML_Chib, file = "T100_ML_Chib.RData")
rm(T100_ML_Chib, T100_obj)

load("T011.RData")
T011_ML_Chib <- MLChib_TVPSSV(Chain = T011_obj, numCores = numCores)
save(T011_ML_Chib, file = "T011_ML_Chib.RData")

load("T101.RData")
T101_ML_Chib <- MLChib_TVPSSV(Chain = T101_obj, numCores = numCores)
save(T101_ML_Chib, file = "T101_ML_Chib.RData")

load("T110.RData")
T110_ML_Chib <- MLChib_TVPSSV(Chain = T110_obj, numCores = numCores)
save(T110_ML_Chib, file = "T110_ML_Chib.RData")

load("T111.RData")
T111_ML_Chib <- MLChib_TVPSSV(Chain = T111_obj, numCores = numCores)
save(T111_ML_Chib, file = "T111_ML_Chib.RData")


setwd("/home/hoanguc3m/Downloads/WP11/ML")
load("G000_ML_Chib.RData"); load("G001_ML_Chib.RData"); load("G010_ML_Chib.RData"); load("G100_ML_Chib.RData")
load("G011_ML_Chib.RData"); load("G101_ML_Chib.RData"); load("G110_ML_Chib.RData"); load("G111_ML_Chib.RData")

load("T000_ML_Chib.RData"); load("T001_ML_Chib.RData"); load("T010_ML_Chib.RData"); load("T100_ML_Chib.RData")
load("T011_ML_Chib.RData"); load("T101_ML_Chib.RData"); load("T110_ML_Chib.RData"); load("T111_ML_Chib.RData")


xtable::xtable(
  round(rbind(c(G000_ML_Chib$CML, G100_ML_Chib$CML, G010_ML_Chib$CML, G001_ML_Chib$CML, G110_ML_Chib$CML, G101_ML_Chib$CML, G011_ML_Chib$CML, G111_ML_Chib$CML),
              c(T000_ML_Chib$CML, T100_ML_Chib$CML, T010_ML_Chib$CML, T001_ML_Chib$CML, T110_ML_Chib$CML, T101_ML_Chib$CML, T011_ML_Chib$CML, T111_ML_Chib$CML)), digits = 2)
)

