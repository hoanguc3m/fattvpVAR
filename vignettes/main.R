library(readxl)
library(Matrix)
library(fattvpVAR)
library(profvis)
library(invgamma)
setwd("/home/hoanguc3m/Downloads/WP11/")
dataraw <- read_excel("/home/hoanguc3m/MEGA/HybridVAR/EconLetter/temp/Data210927.xlsx",
                      col_types = c("text", "numeric", "numeric", "numeric"))
priors <- list(  hyper_ab = 1,
                 hyper_h = 1)

p <- 3 # number of lags
atT <- 821


y <- data.matrix(dataraw[(p+1):atT,c(2:4)])
y0 <- data.matrix(dataraw[1:p,c(2:4)])
#vars::VARselect(y)

# inits <- list(samples = 20000,
#               burnin = 5000,
#               thin = 4,
#               is_tv = c(0,0,0) )

inits <- list(samples = 100000, burnin = 20000, thin = 20)
RhpcBLASctl::blas_set_num_threads(2)

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

setwd("/home/hoanguc3m/Downloads/WP11/")
# library(fatBVARS)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggthemr)
ggthemr('light')

Time <- seq(as.Date("1953/04/01"), as.Date("2021/08/01"), "months")
recessions.df = read.table(textConnection(
  "Peak, Trough
  1953-08-01, 1954-06-01
  1957-09-01, 1958-05-01
  1960-05-01, 1961-03-01
  1969-12-01, 1970-11-01
  1973-11-01, 1975-03-01
  1980-01-01, 1980-07-01
  1981-07-01, 1982-11-01
  1990-07-01, 1991-03-01
  2001-03-01, 2001-11-01
  2007-12-01, 2009-06-01
  2020-03-01, 2020-05-01"), sep=',',
  colClasses=c('Date', 'Date'), header=TRUE)
data_nu <- data.frame(Time = Time, TBILL = dataraw$TBILL3M, S10M3M = dataraw$S10MINUS3M,
                      BAA = dataraw$BAA_SPREAD)

p1 <- ggplot(data_nu, aes(x = Time, y = TBILL)) +
  geom_line(color = "#e3625b") + xlab("TBILL 3M") + ylab("") + ggtitle("Data") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) + theme_bw()

p2 <- ggplot(data_nu, aes(x = Time, y = S10M3M)) +
  geom_line(color = "#e3625b") + xlab("S10-3M") + ylab("") +
  geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) + theme_bw()

p3 <- ggplot(data_nu, aes(x = Time, y = BAA)) +
  geom_line(color = "#e3625b") + xlab("BAA spread") + ylab("") +
  geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5) + theme_bw()

pdf(file='img/data.pdf', width = 9, height = 6)
plot_grid(p1, p2, p3, ncol = 1, align = "v")
#grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()
####################################################################
load("T111.RData")
load("T000.RData")
Nu_mat111 <- T111_obj$store_nu
Nu_mat000 <- T000_obj$store_nu
ndraws <- nrow(Nu_mat000)
varname <- c("TBILL 3M", "S10-3M", "BAA spread")
gg_nu_mat <- function(i){
  data_nu <- data.frame(nu = c(Nu_mat000[,i], Nu_mat111[,i]), Models = c(rep("T000", ndraws), rep("T111", ndraws)))
  max_nu <- max(c(Nu_mat000[,i], Nu_mat111[,i]))
  ggplot(data_nu, aes(x=nu,..density.., fill = Models, colour = Models)) +
    geom_histogram(position = "identity", alpha = 0.8, bins = 30) + xlim(0,max_nu) +
    geom_density(position = "identity", alpha = 0.5) + ylab(expression(nu)) + xlab(varname[i]) + theme_bw() + theme(legend.position="bottom")
}

p1 <- gg_nu_mat(1)
p2 <- gg_nu_mat(2)
p3 <- gg_nu_mat(3)

pdf(file='img/postNu.pdf', width = 12, height = 4)
par(mar=c(2,5,3,1))
grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()
####################################################################
Sigma_h111 <- T111_obj$store_Sigh
Sigma_h000 <- T000_obj$store_Sigh
gg_sigma_mat <- function(i){
  data_nu <- data.frame(nu = c(Sigma_h000[,i], Sigma_h111[,i]), Models = c(rep("T000", ndraws), rep("T111", ndraws)))
  max_nu <- max(c(Sigma_h000[,i], Sigma_h111[,i]))
  ggplot(data_nu, aes(x=nu,..density.., fill = Models, colour = Models)) +
    geom_histogram(position = "identity", alpha = 0.8, bins = 30) + xlim(0,max_nu) +
    geom_density(position = "identity", alpha = 0.5) + ylab(expression(sigma[h]^2)) + xlab(varname[i]) + theme_bw() + theme(legend.position="bottom")
}
p1 <- gg_sigma_mat(1)
p2 <- gg_sigma_mat(2)
p3 <- gg_sigma_mat(3)

pdf(file='img/postSigmah.pdf', width = 12, height = 4)
par(mar=c(2,5,3,1))
grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

####################################################################
Sigma_ab111 <- cbind(T111_obj$store_Sigbeta, T111_obj$store_Sigalp)
varname <- c("b[1]", "B1[1,1]", "B1[1,2]", "B1[1,3]", "B2[1,1]", "B2[1,2]", "B2[1,3]",  "B3[1,1]", "B3[1,2]", "B3[1,3]",
             "b[2]", "B1[2,1]", "B1[2,2]", "B1[2,3]", "B2[2,1]", "B2[2,2]", "B2[2,3]",  "B3[2,1]", "B3[2,2]", "B3[2,3]",
             "b[3]", "B1[3,1]", "B1[3,2]", "B1[3,3]", "B2[3,1]", "B2[3,2]", "B2[3,3]",  "B3[3,1]", "B3[3,2]", "B3[3,3]",
             "A[2,1]", "A[3,1]", "A[3,2]")

gg_sigmaAB_mat <- function(i){
  data_nu <- data.frame(nu = c((Sigma_ab111[,i])) , Models = c(rep("T111", ndraws)))
  max_nu <- (max(c(Sigma_ab111[,i])))
  ggplot(data_nu, aes(x=nu,..density.., fill = Models, colour = Models)) +
    geom_histogram(position = "identity", alpha = 0.8, bins = 30) + xlim(0,max_nu) +
    geom_density(position = "identity", alpha = 0.5) + xlab(bquote(.(varname[i]))) + ylab(expression(sigma^2)) + theme_bw() + theme(legend.position="bottom")
}
l <- list()
for (i in c(1:33)) l[[i]] <- gg_sigmaAB_mat(i)


pdf(file='img/postSigmaAB.pdf', width = 12, height = 3*11)
par(mar=c(2,5,3,1))
plot_grid(l[[1]], l[[11]], l[[21]],
          l[[2]], l[[12]], l[[22]],
          l[[3]], l[[13]], l[[23]],
          l[[4]], l[[14]], l[[24]],
          l[[5]], l[[15]], l[[25]],
          l[[6]], l[[16]], l[[26]],
          l[[7]], l[[17]], l[[27]],
          l[[8]], l[[18]], l[[28]],
          l[[9]], l[[19]], l[[29]],
          l[[10]], l[[20]], l[[30]],
          l[[31]], l[[32]], l[[33]],
          ncol = 3, byrow = TRUE)
dev.off()

####################################################################
ab111_mean <- cbind( apply(T111_obj$store_beta, MARGIN = c(2,3), FUN = mean),
                     apply(T111_obj$store_alp, MARGIN = c(2,3), FUN = mean))
ab111_q10 <- cbind( apply(T111_obj$store_beta, MARGIN = c(2,3), FUN = quantile, probs = c(0.1)),
                     apply(T111_obj$store_alp, MARGIN = c(2,3), FUN = quantile, probs = c(0.1)))
ab111_q90 <- cbind( apply(T111_obj$store_beta, MARGIN = c(2,3), FUN = quantile, probs = c(0.9)),
                    apply(T111_obj$store_alp, MARGIN = c(2,3), FUN = quantile, probs = c(0.9)))

gg_TVPAB_mat <- function(i){
  Time <- tail(seq(as.Date("1953/04/01"), as.Date("2021/08/01"), "months"), 818)

  data_nu <- data.frame(Time = Time, AB_mean = ab111_mean[,i])

  data_nu1 <- data.frame(Time = c(Time, rev(Time)) ,
                         AB_UL = c(ab111_q10[,i], rev(ab111_q90[,i])))
  ggplot() + geom_line(data=data_nu, mapping=aes(x=Time, y=AB_mean), color = "#ff0c00") +
    geom_polygon(data=data_nu1, mapping=aes(x=Time, y=AB_UL), fill = "#5ba6d6", alpha = 0.5) +
    xlab(bquote(.(varname[i]))) + ylab("") + theme_bw() + theme(legend.position="bottom") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_rect(data=recessions.df, inherit.aes=F, aes(xmin=Peak, xmax=Trough, ymin=-Inf, ymax=+Inf), fill='darkgray', alpha=0.5)

}

l <- list()
for (i in c(1:33)) l[[i]] <- gg_TVPAB_mat(i)


pdf(file='img/postAB.pdf', width = 12, height = 3*11)
par(mar=c(2,5,3,1))
plot_grid(l[[1]], l[[11]], l[[21]],
          l[[2]], l[[12]], l[[22]],
          l[[3]], l[[13]], l[[23]],
          l[[4]], l[[14]], l[[24]],
          l[[5]], l[[15]], l[[25]],
          l[[6]], l[[16]], l[[26]],
          l[[7]], l[[17]], l[[27]],
          l[[8]], l[[18]], l[[28]],
          l[[9]], l[[19]], l[[29]],
          l[[10]], l[[20]], l[[30]],
          l[[31]], l[[32]], l[[33]],
          ncol = 3, byrow = TRUE)
dev.off()
