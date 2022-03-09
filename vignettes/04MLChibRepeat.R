library(fattvpVAR)
library(RhpcBLASctl)
numCores = 4
blas_set_num_threads(numCores)
setwd("/Backup/Ongoing/WP11/")
setwd("/home/hoanguc3m/Downloads/WP11")
load("~/Dropbox/WP11/Code/fattvpVAR/data/Spread.RData")

set.seed(NULL)



##########################################################################
# inits$is_tv = c(0,0,0); G000_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G000_ML_Chib <- MLChib_TVPGSV(Chain = G000_obj, numCores = numCores)
# save(G000_ML_Chib, file = paste("G000_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G000_ML_Chib, G000_obj)


# inits$is_tv = c(0,0,1); G001_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G001_ML_Chib <- MLChib_TVPGSV(Chain = G001_obj, numCores = numCores)
# save(G001_ML_Chib, file = paste("G001_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G001_ML_Chib, G001_obj)

# inits$is_tv = c(0,1,0); G010_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G010_ML_Chib <- MLChib_TVPGSV(Chain = G010_obj, numCores = numCores)
# save(G010_ML_Chib, file = paste("G010_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G010_ML_Chib, G010_obj)
#
# inits$is_tv = c(1,0,0); G100_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G100_ML_Chib <- MLChib_TVPGSV(Chain = G100_obj, numCores = numCores)
# save(G100_ML_Chib, file = paste("G100_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G100_ML_Chib, G100_obj)
#
# inits$is_tv = c(0,1,1); G011_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G011_ML_Chib <- MLChib_TVPGSV(Chain = G011_obj, numCores = numCores)
# save(G011_ML_Chib, file = paste("G011_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G011_ML_Chib, G011_obj)
#
# inits$is_tv = c(1,0,1); G101_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G101_ML_Chib <- MLChib_TVPGSV(Chain = G101_obj, numCores = numCores)
# save(G101_ML_Chib, file = paste("G101_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G101_ML_Chib, G101_obj)
#
# inits$is_tv = c(1,1,0); G110_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G110_ML_Chib <- MLChib_TVPGSV(Chain = G110_obj, numCores = numCores)
# save(G110_ML_Chib, file = paste("G110_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G110_ML_Chib, G110_obj)
#
#
#
# inits$is_tv = c(1,1,1); G111_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
# G111_ML_Chib <- MLChib_TVPGSV(Chain = G111_obj, numCores = numCores)
# save(G111_ML_Chib, file = paste("G111_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(G111_ML_Chib, G111_obj)
#
# ##########################################################################
#
# inits$is_tv = c(0,0,0); T000_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
# T000_ML_Chib <- MLChib_TVPSSV(Chain = T000_obj, numCores = numCores)
# save(T000_ML_Chib, file = paste("T000_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(T000_ML_Chib, T000_obj)
#
#
# inits$is_tv = c(0,0,1); T001_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
# T001_ML_Chib <- MLChib_TVPSSV(Chain = T001_obj, numCores = numCores)
# save(T001_ML_Chib, file = paste("T001_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(T001_ML_Chib, T001_obj)
#
# inits$is_tv = c(0,1,0); T010_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
# T010_ML_Chib <- MLChib_TVPSSV(Chain = T010_obj, numCores = numCores)
# save(T010_ML_Chib, file = paste("T010_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(T010_ML_Chib, T010_obj)
#
# inits$is_tv = c(1,0,0); T100_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
# T100_ML_Chib <- MLChib_TVPSSV(Chain = T100_obj, numCores = numCores)
# save(T100_ML_Chib, file = paste("T100_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(T100_ML_Chib, T100_obj)
#
# inits$is_tv = c(0,1,1); T011_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
# T011_ML_Chib <- MLChib_TVPSSV(Chain = T011_obj, numCores = numCores)
# save(T011_ML_Chib, file = paste("T011_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(T011_ML_Chib, T011_obj)
#
# inits$is_tv = c(1,0,1); T101_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
# T101_ML_Chib <- MLChib_TVPSSV(Chain = T101_obj, numCores = numCores)
# save(T101_ML_Chib, file = paste("T101_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(T101_ML_Chib, T101_obj)
#
# inits$is_tv = c(1,1,0); T110_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
# T110_ML_Chib <- MLChib_TVPSSV(Chain = T110_obj, numCores = numCores)
# save(T110_ML_Chib, file = paste("T110_ML_Chib_", Sys.time() , ".RData", sep = ""))
# rm(T110_ML_Chib, T110_obj)


inits$is_tv = c(1,1,1); T111_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
T111_ML_Chib <- MLChib_TVPSSV(Chain = T111_obj, numCores = numCores)
save(T111_ML_Chib, file = paste("T111_ML_Chib_", Sys.time() , ".RData", sep = ""))
rm(T111_ML_Chib, T111_obj)

#
#
# setwd("/home/hoanguc3m/Downloads/WP11/ML")
# load("G000_ML_Chib.RData"); load("G001_ML_Chib.RData"); load("G010_ML_Chib.RData"); load("G100_ML_Chib.RData")
# load("G011_ML_Chib.RData"); load("G101_ML_Chib.RData"); load("G110_ML_Chib.RData"); load("G111_ML_Chib.RData")
#
# load("T000_ML_Chib.RData"); load("T001_ML_Chib.RData"); load("T010_ML_Chib.RData"); load("T100_ML_Chib.RData")
# load("T011_ML_Chib.RData"); load("T101_ML_Chib.RData"); load("T110_ML_Chib.RData"); load("T111_ML_Chib.RData")
#
#
# xtable::xtable(
#   round(rbind(c(G000_ML_Chib$CML, G100_ML_Chib$CML, G010_ML_Chib$CML, G001_ML_Chib$CML, G110_ML_Chib$CML, G101_ML_Chib$CML, G011_ML_Chib$CML, G111_ML_Chib$CML),
#               c(T000_ML_Chib$CML, T100_ML_Chib$CML, T010_ML_Chib$CML, T001_ML_Chib$CML, T110_ML_Chib$CML, T101_ML_Chib$CML, T011_ML_Chib$CML, T111_ML_Chib$CML)), digits = 2)
# )

