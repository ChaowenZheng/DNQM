rm(list=ls())

# install and load packages  "xlsx", 
libraries = c("latex2exp", "xtable", "quantreg", "lubridate", "bigmemory", "igraph", "WRTDStidal", "Matrix")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

folder = "C:/Users/sprin/Desktop/Application20211108/github"
setwd(folder) 
getwd()

source("netestimator.R")

yeardata = 2016
taus = seq(0.1, 0.9, 0.1)

##################
YNetData = readRDS("YNetData")
Ymat = YNetData$Ymat
Z = YNetData$Z
Fat = YNetData$sFat


Nfirm = nrow(Ymat)   # firm 943
Time = nrow(Fat)
### network 

WCity = YNetData$WCity
W = WCity   #  Network by the same city
Wdensity = nnzero(W)/(Nfirm^2-Nfirm)
Wdensity



coef = matrix(0, nrow = 4+3+4+1, ncol = length(taus))
se = matrix(0, nrow = 4+3+4+1, ncol = length(taus))
t_val = matrix(0, nrow = 4+3+4+1, ncol = length(taus))
p_val = matrix(0, nrow = 4+3+4+1, ncol = length(taus))

Lagcoef = matrix(0, nrow = 4+3+4, ncol = length(taus))
Lagse = matrix(0, nrow = 4+3+4, ncol = length(taus))
Lagt_val = matrix(0, nrow = 4+3+4, ncol = length(taus))
Lagp_val = matrix(0, nrow = 4+3+4, ncol = length(taus))

Fcoef = matrix(0, nrow = 4+3+1, ncol = length(taus))
Fse = matrix(0, nrow = 4+3+1, ncol = length(taus))
Ft_val = matrix(0, nrow = 4+3+1, ncol = length(taus))
Fp_val = matrix(0, nrow = 4+3+1, ncol = length(taus))

FLagcoef = matrix(0, nrow = 4+3, ncol = length(taus))
FLagse = matrix(0, nrow = 4+3, ncol = length(taus))
FLagt_val = matrix(0, nrow = 4+3, ncol = length(taus))
FLagp_val = matrix(0, nrow = 4+3, ncol = length(taus))

Netallout = list()
FNetallout = list()
Lagnetallout = list()
FLagnetallout = list()

goodf = matrix(0, nrow = 5, ncol = length(taus))

for (q in 1:length(taus))
{
  # q = 1
  tau = taus[q]
  cat(" W = W1,  tau ", tau, "\n")
  
  # IV
  Ymat1 = Ymat[, -1]
  M21 = W%*%W%*%Ymat1 
  M31 = W%*%W%*%W%*%Ymat1
  M = cbind(as.vector(M21), as.vector(M31))
  
  Netqrout = netQuantile(Ymat, W, Z, M, tau, Fat)
  Netallout[[q]] = Netqrout
  coef[, q] = Netqrout$b[, 1]
  se[, q]   = Netqrout$b[, 2]
  t_val[, q] = Netqrout$b[, 3]
  p_val[, q] = Netqrout$b[, 4]
  
#  non simutaneous elememt NQAR 
  Lagnetout = lagnetQuantile(Ymat, W, Z, tau, Fat)  
  Lagnetallout[[q]] = Lagnetout
  Lagcoef[, q] = Lagnetout$coefficients[, 1]
  Lagse[, q]   = Lagnetout$coefficients[, 2]
  Lagt_val[, q] = Lagnetout$coefficients[, 3]
  Lagp_val[, q] = Lagnetout$coefficients[, 4]
  
  # no common factor
  FNetqrout = netQuantile(Ymat, W, Z, M, tau)
  FNetallout[[q]] = FNetqrout
  Fcoef[, q] = FNetqrout$b[, 1]
  Fse[, q]   = FNetqrout$b[, 2]
  Ft_val[, q] = FNetqrout$b[, 3]
  Fp_val[, q] = FNetqrout$b[, 4]
  
  #  non simutaneous elememt no common factor
  FLagnetout = lagnetQuantile(Ymat, W, Z, tau)  
  FLagnetallout[[q]] = FLagnetout
  FLagcoef[, q] = FLagnetout$coefficients[, 1]
  FLagse[, q]   = FLagnetout$coefficients[, 2]
  FLagt_val[, q] = FLagnetout$coefficients[, 3]
  FLagp_val[, q] = FLagnetout$coefficients[, 4]
  
  # goodness of fit
  goodf[1, q] = goodfit(Netqrout$res, Lagnetout$res, tau)
  goodf[2, q] = goodfit(Netqrout$res, FNetqrout$res, tau)
  goodf[3, q] = goodfit(Netqrout$res, FLagnetout$res, tau)
  goodf[4, q] = goodfit(Lagnetout$res, FLagnetout$res, tau)
  goodf[5, q] = goodfit(FNetqrout$res, FLagnetout$res, tau)
  
}


figname = paste(folder, '/Figure/', yeardata, 'figparWCityallm', 'N', Nfirm, '.pdf', sep = "")
pdf(figname, paper="a4", height = 0, width = 0)
Plotnetivrq(Netallout, parm = seq(2, 12, 1))
dev.off() 






