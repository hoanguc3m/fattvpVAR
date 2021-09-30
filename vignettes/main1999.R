library(readxl)
library(Matrix)
library(fattvpVAR)
library(profvis)
library(invgamma)

dataraw <- read_excel("/home/ORUNET.ORU.SE/hgnn/Downloads/TVPVAR/tmp/Data210927.xlsx",
                      col_types = c("text", "numeric", "numeric", "numeric"))
priors <- list(  hyper_ab = 1,
                 hyper_h = 1)

p <- 3 # number of lags
atT <- 561
dataraw$Time[atT]
setwd("/home/ORUNET.ORU.SE/hgnn/Downloads/TVPVAR/1999")


y <- data.matrix(dataraw[(p+1):atT,c(2:4)])
y0 <- data.matrix(dataraw[1:p,c(2:4)])
#vars::VARselect(y)

# inits <- list(samples = 20000,
#               burnin = 5000,
#               thin = 4,
#               is_tv = c(0,0,0) )

inits <- list(samples = 100000, burnin = 20000, thin = 20)
RhpcBLASctl::blas_set_num_threads(3)

####################################################################
{
inits$is_tv = c(0,0,0); G000_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G000_obj, file = "G000.RData")


inits$is_tv = c(0,0,1); G001_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G001_obj, file = "G001.RData")

inits$is_tv = c(0,1,0); G010_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G010_obj, file = "G010.RData")

inits$is_tv = c(1,0,0); G100_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G100_obj, file = "G100.RData")

inits$is_tv = c(1,1,0); G110_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G110_obj, file = "G110.RData")

inits$is_tv = c(0,1,1); G011_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G011_obj, file = "G011.RData")

inits$is_tv = c(1,0,1); G101_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G101_obj, file = "G101.RData")

inits$is_tv = c(1,1,1); G111_obj <- GaussTVPSV(y, y0, p, priors, inits)
save(G111_obj, file = "G111.RData")
}
####################################################################
{
inits$is_tv = c(0,0,0); T000_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T000_obj, file = "T000.RData")

inits$is_tv = c(0,0,1); T001_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T001_obj, file = "T001.RData")

inits$is_tv = c(0,1,0); T010_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T010_obj, file = "T010.RData")

inits$is_tv = c(1,0,0); T100_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T100_obj, file = "T100.RData")

inits$is_tv = c(1,1,0); T110_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T110_obj, file = "T110.RData")

inits$is_tv = c(0,1,1); T011_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T011_obj, file = "T011.RData")

inits$is_tv = c(1,0,1); T101_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T101_obj, file = "T101.RData")

inits$is_tv = c(1,1,1); T111_obj <- StudentTVPSV(y, y0, p, priors, inits)
save(T111_obj, file = "T111.RData")
}
####################################################################
