library(readxl)
library(Matrix)
library(fattvpVAR)
library(profvis)
library(invgamma)
dataraw <- read_excel("/home/hoanguc3m/MEGA/HybridVAR/EconLetter/temp/Data210324.xlsx",
                      col_types = c("text", "numeric", "numeric", "numeric"))
inits <- list(samples = 10000,
              burnin = 1000,
              thin = 1,
              is_tv = c(1,1,1) )
priors <- list(  hyper_ab = 1,
                 hyper_h = 1)

p <- 3 # number of lags
atT <- 815

y <- data.matrix(dataraw[(p+1):atT,c(2:4)])
y0 <- data.matrix(dataraw[1:p,c(2:4)])
G111_obj <- GaussTVPSV(y, y0, p, priors, inits)
T111_obj <- StudentTVPSV(y, y0, p, priors, inits)

