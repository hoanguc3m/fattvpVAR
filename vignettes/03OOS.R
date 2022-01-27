library(readxl)
library(Matrix)
library(fattvpVAR)
library(profvis)
library(invgamma)
setwd("/home/hoanguc3m/Downloads/WP11/")
dataraw <- read_excel("/home/hoanguc3m/MEGA/HybridVAR/EconLetter/temp/Data210927.xlsx",
                      col_types = c("text", "numeric", "numeric", "numeric"))
p <- 3 # number of lags
y_raw = dataraw[,c(2:4)]
c(dataraw[561,1], dataraw[797,1])


RhpcBLASctl::blas_set_num_threads(1)
numCores <- 16

#####################################################
is_tv = c(0,0,0)
recursive000 <- parallel::mclapply(561:797,
                              FUN = function(t_start) {
                                recursive_seperate(y_raw, is_tv = is_tv, t_start = t_start,
                                                   t_pred = 24, K = ncol(y_raw), p = p,
                                                   dist = "Gaussian", outname = NULL)
                              }, mc.cores = numCores)
save(recursive000, "RecursiveG000.RData")

#####################################################
is_tv = c(0,0,1)
recursive001 <- parallel::mclapply(561:797,
                                   FUN = function(t_start) {
                                     recursive_seperate(y_raw, is_tv = is_tv, t_start = t_start,
                                                        t_pred = 24, K = ncol(y_raw), p = p,
                                                        dist = "Gaussian", outname = NULL)
                                   }, mc.cores = numCores)
save(recursive001, "RecursiveG001.RData")

#####################################################
is_tv = c(1,1,1)
recursive111 <- parallel::mclapply(561:797,
                                   FUN = function(t_start) {
                                     recursive_seperate(y_raw, is_tv = is_tv, t_start = t_start,
                                                        t_pred = 24, K = ncol(y_raw), p = p,
                                                        dist = "Gaussian", outname = NULL)
                                   }, mc.cores = numCores)
save(recursive111, "RecursiveG111.RData")


#####################################################


#####################################################
is_tv = c(0,0,0)
recursive000 <- parallel::mclapply(561:797,
                                   FUN = function(t_start) {
                                     recursive_seperate(y_raw, is_tv = is_tv, t_start = t_start,
                                                        t_pred = 24, K = ncol(y_raw), p = p,
                                                        dist = "Student", outname = NULL)
                                   }, mc.cores = numCores)
save(recursive000, "RecursiveT000.RData")

#####################################################
is_tv = c(0,0,1)
recursive001 <- parallel::mclapply(561:797,
                                   FUN = function(t_start) {
                                     recursive_seperate(y_raw, is_tv = is_tv, t_start = t_start,
                                                        t_pred = 24, K = ncol(y_raw), p = p,
                                                        dist = "Student", outname = NULL)
                                   }, mc.cores = numCores)
save(recursive001, "RecursiveT001.RData")

#####################################################
is_tv = c(1,1,1)
recursive111 <- parallel::mclapply(561:797,
                                   FUN = function(t_start) {
                                     recursive_seperate(y_raw, is_tv = is_tv, t_start = t_start,
                                                        t_pred = 24, K = ncol(y_raw), p = p,
                                                        dist = "Student", outname = NULL)
                                   }, mc.cores = numCores)
save(recursive111, "RecursiveT111.RData")


#####################################################
