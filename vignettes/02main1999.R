library(readxl)
library(Matrix)
library(fattvpVAR)
library(profvis)
library(invgamma)
# setwd("/home/hoanguc3m/Downloads/WP11/")
# dataraw <- read_excel("/home/hoanguc3m/MEGA/HybridVAR/EconLetter/temp/Data210927.xlsx",
#                       col_types = c("text", "numeric", "numeric", "numeric"))
# priors <- list(  hyper_ab = 1,
#                  hyper_h = 1)
#
# p <- 3 # number of lags
# atT <- 561
# dataraw$Time[atT]
#
#
# y <- data.matrix(dataraw[(p+1):atT,c(2:4)])
# y0 <- data.matrix(dataraw[1:p,c(2:4)])
# priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE)
# inits <- list(samples = 100000, burnin = 20000, thin = 20)
load("~/Dropbox/WP11/Code/fattvpVAR/data/Spread.RData")
atT <- 561
dataraw$Time[atT]
y <- data.matrix(dataraw[(p+1):atT,c(2:4)])
y0 <- data.matrix(dataraw[1:p,c(2:4)])
priors <- get_prior_minnesota(y = y, p = p, intercept=TRUE)
inits <- list(samples = 100000, burnin = 20000, thin = 20)


####################################################################
{
inits$is_tv = c(0,0,0); G000_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G000_obj, file = "G000.RData")


inits$is_tv = c(0,0,1); G001_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G001_obj, file = "G001.RData")

inits$is_tv = c(0,1,0); G010_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G010_obj, file = "G010.RData")

inits$is_tv = c(1,0,0); G100_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G100_obj, file = "G100.RData")

inits$is_tv = c(1,1,0); G110_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G110_obj, file = "G110.RData")

inits$is_tv = c(0,1,1); G011_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G011_obj, file = "G011.RData")

inits$is_tv = c(1,0,1); G101_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G101_obj, file = "G101.RData")

inits$is_tv = c(1,1,1); G111_obj <- fitTVPGaussSV(y, y0, p, priors, inits)
save(G111_obj, file = "G111.RData")
}
####################################################################
{
inits$is_tv = c(0,0,0); T000_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T000_obj, file = "T000.RData")

inits$is_tv = c(0,0,1); T001_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T001_obj, file = "T001.RData")

inits$is_tv = c(0,1,0); T010_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T010_obj, file = "T010.RData")

inits$is_tv = c(1,0,0); T100_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T100_obj, file = "T100.RData")

inits$is_tv = c(1,1,0); T110_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T110_obj, file = "T110.RData")

inits$is_tv = c(0,1,1); T011_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T011_obj, file = "T011.RData")

inits$is_tv = c(1,0,1); T101_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T101_obj, file = "T101.RData")

inits$is_tv = c(1,1,1); T111_obj <- fitTVPStudentSV(y, y0, p, priors, inits)
save(T111_obj, file = "T111.RData")
}
####################################################################
