rm(list=ls())
#################################################################################################################
#########################                  Read me                 --------------------------------------------
#################################################################################################################
#We first define many functions, and you could read start from the MAIN CODE at around line  550
#The code before 550 are mainly about define functions for generating data and estimation
#The code after line 847 are mainly about storing the results and making latex tables
#After changing the working directory, you should able to run it directly 
getwd()
#change to your own working direct
setwd("")
# install and load packages

libraries = c("MASS", "Matrix", "poweRlaw", "methods","xtable", "quantreg")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)



#################################################################################################################
#########################                  Functions                 --------------------------------------------
#################################################################################################################

# 1.######################################Estimation Function ---------------------------------------------------
##### IV directly
#The input:
#Ymat:N *T Dependent variable
#BN: The burn-in sample for generating dynamic data; In empirical anlysis, it could be set as 1
#W:spatial weighting matrix
##X:N *T Inependent variable
#Z:time invariant variable
#Fatt:factors
#tau: quantile level
#SR:range of the spatial coefficient when doing grid search
#IVs: Specify the IVs to be used.


netQuantileIVdir<-function(Ymat,BN, W, Z, Fatt, tau, SR,IVs) 
{
  WYmat = W%*%Ymat[,(BN+1):ncol(Ymat)]
  Time = ncol(WYmat)
  #Wide to long: spatial variable
  D<-as.vector(WYmat)
  #Dynamic
  Ymatlag<-Ymat[,BN:(ncol(Ymat)-1)]
  dim(Ymatlag)
  #spatial lag
  WYmatlag = W%*%Ymatlag
  dim(WYmatlag)
  
  #Factors
  Fat = Fatt[-1, ]
  f = cbind(Fat[(BN+1):nrow(Fat), 1], Fat[(BN):(nrow(Fat)-1), 1],
            Fat[(BN+1):nrow(Fat), 2], Fat[(BN):(nrow(Fat)-1), 2])
  Fatlist = lapply(1:4, function(i) {Ft = t(f[, i]%x%matrix(1, 1, nrow(W)))
  return(as.vector(Ft)) })
  
  # all regressors
  X_all = cbind(as.vector(Ymatlag),   # gamma2 lagged network effects
                as.vector(WYmatlag),   # gamma3 lagged y
                do.call("rbind", rep(list(Z), Time)), # alpha1-alpha5  5
                do.call("cbind", Fatlist)) # beta1-beta4   4  
  
  ## IV WX and W^2X
  if(IVs=="wx1"){
    IV1 = W%*%X
    M = as.matrix(as.vector(IV1))
  } 
  if(IVs=="wx12"){
    IV1 = W%*%X
    IV2 = W%*%W%*%X
    M = cbind(as.vector(IV1),as.vector(IV2))
    dim(M)
  } 
  if(IVs=="wy2"){
    IV1 = W%*%WYmatlag
    M = as.matrix(as.vector(IV1))
  }
  if(IVs=="wy3"){
    IV1 = W%*%W%*%WYmatlag
    M = as.matrix(as.vector(IV1))
  }
  if(IVs=="wy23"){
    IV1 = W%*%WYmatlag
    IV2 = W%*%W%*%WYmatlag
    M = cbind(as.vector(IV1),as.vector(IV2))
    dim(M)
  }
  
  
  # M is IV, XIV including X and IV
  XIV = cbind(M,X_all)
  dim(XIV)
  # dependent variable
  Yvec = as.vector(Ymat[,(BN+1):ncol(Ymat)])
  
  
  ## Step 1 set the range of SR or alpha parameter   # tail(alpha, 3)
 
  if(is.null(SR)) {
    #this method could be problematic sometimes, depending on the lengths of each step
    ## the initial value 
    #  independent variable X without IV
    XIVnon = cbind(D, X_all)           
    
    qrnonIV = rq(Yvec ~ XIVnon, tau)
    qrnonIV_coef = summary.rq(qrnonIV, se = "boot", hs = F)$coefficients 
    gamma1_IVnon = qrnonIV_coef[2, 1]
    std_IVnon = qrnonIV_coef[2, 2]
    alpha = seq(gamma1_IVnon - 4*std_IVnon, 
                gamma1_IVnon + 4*std_IVnon, 0.1*std_IVnon)
    
  } else {
    alpha = SR
  }
  
  
  ## Step 2 get alpha by the minimum value of IV
  gnorm = lapply(alpha, function(t)
  {
    # cat(" ", t)
    dat   = data.frame(Y = Yvec - t*as.vector(WYmat), XIV)
    resrq = rq(Y~., tau, data = dat)
    infer_coef = resrq$coefficients   # boot
    g = norm(infer_coef[2:(1+ncol(M))])
    return(g)
  })
  gnorm  = do.call("rbind", gnorm)
  gamma1 = alpha[which(gnorm == min(gnorm))]
  #cat("Min: ", gamma1,"Quantile:",tau)
  
  ## Step 3 plug in alpha and get other estimators
  dat   = data.frame(Y = Yvec - gamma1*D,X_all)
  resrq = rq(Y~., tau, data = dat)
  
  
  ## Compute the s.e. 
  n = nrow(dat)

  
  #Calculating the density function
  se<-matrix(0,ncol(X_all)+2,4)
  
  #Bandwidth
  bh_ind<-0
  Xmat = cbind(1, M,X_all)
  Xmatt = cbind(1, D, X_all)

  
  for(Band_type in c("HS","HB")){
    bh_ind<-bh_ind+1
    if (Band_type=="HS"){#Hall and Sheather(1988)
      scale_range<-c(1,3) #Following Koenekr and Xiao, 2006. we may consider scaled bandwidth
    } else {#Bofinger(1975)
      scale_range<-c(1,0.6) #Following Koenekr and Xiao, 2006. we may consider scaled bandwidth
    }
    sh_ind<-0
    for (scale_h in scale_range){
      sh_ind<-sh_ind+1
      if (Band_type=="HS"){
        #First, the so called difference quotient method
        hf<-n^(-1/3)*(qnorm(0.975))^(2/3)*(1.5*(dnorm(qnorm(tau))^2)/(2*(qnorm(tau))^2+1))^(1/3)
      } else {
        hf<-n^(-1/5)*(4.5*(dnorm(qnorm(tau))^4)/((2*(qnorm(tau))^2+1)^2))^(1/5)
      }
      
      hf<-scale_h*hf
      #print(hf)
      
      if(tau+hf>=1 |tau-hf<=0){
        hf<-min(abs(tau-0.0001),abs(1-tau-0.0001)) #0.0001 is arbitry selected to ensure the quantile is well defined (within (0,1))
        den<-2*hf/(quantile(resrq$residuals,tau+hf)-quantile(resrq$residuals,tau-hf))
      } else {
        den<-2*hf/(quantile(resrq$residuals,tau+hf)-quantile(resrq$residuals,tau-hf))
      }
      
      Ome1 = t(Xmatt)%*%(Xmat)*den
      Ome0 = crossprod(Xmat)
      covQ = solve(Ome1%*%solve(Ome0)%*%t(Ome1))*tau*(1-tau)
      se[,(bh_ind-1)*2+sh_ind]<-sqrt(abs(diag(covQ)))
    }
  }
  
  ##########Estimation without using IV
  
  XIVnon = cbind(D, X_all)           
  
  qrnonIV = rq(Yvec ~ XIVnon, tau)
  estIVnon<-summary.rq(qrnonIV,se="iid")$coefficients[,1] 
  seIVnon<-summary.rq(qrnonIV,se="iid")$coefficients[,2] 
  
  ## output
  return(list(coefficients = c(resrq$coefficients[1], gamma1, resrq$coefficients[2:(ncol(X_all)+1)]), 
              se = se,estIVno=estIVnon,seIVno=seIVnon))
  
}




# 2. #####################################Data Simulation Functions -------------------------------------------------

#Functions for generating parameters


fphi<-function(u)
{
  phi = 0.2+0.1*u
  return(phi)
}

fphi1<-function(u)
{
  phi1 = 0.4*exp(u)/(1+exp(u))
  return(phi1)
}
fphi2<-function(u)
{
  phi2 = 0.2+0.1*u
  return(phi2)
}


get.SR<-function(u)
{
  SR = seq(phi1(u) - 0.05, phi1(u) + 0.05, 0.01)
  return(SR)
}

##for time invariant and individual invariant variable

get.gamma.nodal<-function(u)
{
  gam = c(0.5*pnorm(u),
          0.3*pgamma(u, shape = 1, scale = 2),
          0.2*pgamma(u, shape = 2, scale = 2),
          0.25*pgamma(u, shape = 3, scale = 2),
          0.2*pgamma(u, shape = 2, scale = 1))
  return(gam)
}


get.nodal<-function(u, Z0)
{#matrix  (a1,a2,...aN)'\times vector (b1,b2,...bN)=(B1a1,B2a2,...bNaN)'
  gList = list()
  gList[[1]] = 0.5*pnorm(u)*Z[, 1]
  gList[[2]] = 0.3*pgamma(u, shape = 1, scale = 2)*Z[, 2]
  gList[[3]] = 0.2*pgamma(u, shape = 2, scale = 2)*Z[, 3]
  gList[[4]] = 0.25*pgamma(u, shape = 3, scale = 2)*Z[, 4]
  gList[[5]] = 0.2*pgamma(u, shape = 2, scale = 1)*Z[, 5]
  nodal = Reduce("+", gList) 
  return(nodal)
}


get.gamma.exog<-function(u)
{
  gam = c(0.1*pnorm(u),
          0.2*pnorm(u),
          0.1*pgamma(u, shape = 2, scale = 2),
          0.3*pgamma(u, shape = 2, scale = 2))
  return(gam)
}



get.exog<-function(u, Fat, FLag = 1)
{ 
  #Sine F is time invariant (the same for every individual at the same time)
  f = cbind(Fat[2:nrow(Fat), 1], Fat[1:(nrow(Fat)-1), 1],
            Fat[2:nrow(Fat), 2], Fat[1:(nrow(Fat)-1), 2])
  Fatlist = lapply(1:(2*(FLag + 1)), function(i) t(f[, i]%x%matrix(1, 1, nrow(u))))
  fList = list()
  fList[[1]] = 0.1*pnorm(u)*Fatlist[[1]]
  fList[[2]] = 0.2*pnorm(u)*Fatlist[[2]]
  fList[[3]] = 0.1*pgamma(u, shape = 2, scale = 2)*Fatlist[[3]]
  fList[[4]] = 0.3*pgamma(u, shape = 2, scale = 2)*Fatlist[[4]]
  exog = Reduce("+", fList)
  return(exog)
}

#Functions for generating data

#simu.data<-function(W, Time ,X, Z,Fat, type)
simu.data<-function(W, Time, Z,Fat, type)
{
  W<-as.matrix(W)
  U = matrix(runif(nrow(W)*Time, 0, 1), nrow = nrow(W))
  if (type == "Z")
  {
    epsilon = qnorm(U)
  }
  if (type == "T")
  {
    epsilon = qt(U,5)
  }
  if (type == "chi")
  {
    epsilon = qchisq(U,3)
  }
  
  #Spatial Coef
  phi = fphi(U)
  mean(phi)
  phi1 = fphi1(U)
  phi2 = fphi2(U)

  #Factors:F
  exog = get.exog(epsilon, Fat)
  
  #Individual invariant:Z
  nodal = get.nodal(epsilon, Z)
 
  #Generating the data
  Ymat<-matrix(0,N,Time)
  for (i in 2:Time){
    Ymat[,i] = as.matrix(solve(diag(1, nrow(W), nrow(W)) - diag(phi[, i])%*%W)%*%(epsilon[,i]+diag(phi1[, i])%*%Ymat[,i-1]+diag(phi2[, i])%*%W%*%Ymat[,i-1]+exog[,i]+nodal[,i]))
  }
  return(Ymat)
}

#3.################################## For generating Spatial weighting matrix ----------------------------------

# simulate the simple asymatic matrix: for each individual, there are around 2h neighbours
getAsymW<-function(N,h)
{
  W<-matrix(0,N,N)
  for (i in 1:N){
    if((i+h) <= N){
      W[i,i:(i+h)]<-1
    }
    else
    {
      #W[i,i]<-1
      W[i,i:ncol(W)]<-1
    }
  }
  W[lower.tri(W)] = t(W)[lower.tri(W)]
  diag(W)<-0
  rs<-matrix(rowSums(W),N,1)
  for (i in 1:N){
    W[i,]<-W[i,]/rs[i,]
  }
  return(W)
}



getDyadW<-function(N, delta, normalize = T)
{
  A = matrix(0, nrow = N, ncol = N)
  
  ############### mutual follower
  #Algorithm: First randomly pick N*delta elements to be assigned 1
  #save their index; and then assign their symmetric elements to be 1 (mutual follower)
  ind = which(upper.tri(A), arr.ind = T)
  #Randomly pick N*delta elements
  indM = ind[sample(1:nrow(ind), N*delta), ]
  A[indM] = 1
  A[indM[, 2:1]] = 1
  
  ############### single relationship
  #Algorithm: First randomly pick N^delta elements to be assigned 1
  #save their index; and then split them evenly to 2, change the order of the index of one group
  indl = which(A==0&upper.tri(A), arr.ind = T)
  indS = indl[sample(1:nrow(indl), N^1.2), ]
  tmp = sample(1:nrow(indS), floor(N^1.2/2))
  indS[tmp, ] = indS[tmp, 2:1]
  A[indS] = 1
  diag(A) = 0
  
  #In case some rowsum is zero
  ind = which(rowSums(A)==0)                                                                          
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 2)] = 1                           
  }
  
  if(!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}


getBlockW<-function(N, Nblock, normalize = T, delta = 0.3)                                                
{
  #Algorithm for generating:
  #1. Define the block and the non-block elements;
  #2.mList is for element within the block
  #3.offDiag are for elements in the between block
  if (N%%Nblock==0){                                                                                  
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)                   
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.3*N^(-delta)),                     
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)        
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    #The Binomial Distribution: size: number of trials;prob:probability of success on each trial.
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.3*N^{-delta}),                
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.3*N^{-delta}),              
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)                                                                       
  offDiag = which(isDiag == 0, arr.ind = T)                                                           
  #If the block has an island that is not connected with anybody else, we then try to avoid this.
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                                       
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  #
  ind = which(rowSums(bA)==0)                                                                         
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                                                             
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                                                                
  W = as(W, "dgCMatrix")
  return(as.matrix(W))
}



getPowerLawW<-function(N, alpha, normalize = T)                                                      
{
  #The probability of taking large values decreases sharply
  Nfollowers = rpldis(N, 1, alpha)                                                                    
  A = sapply(Nfollowers, function(n) {                                                                
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                          
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 4)] = 1                                                              
  }
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  W = as(W, "dgCMatrix")
  return(W)
}

norm<-function(x){sqrt(sum(x^2))}


####################################################################################################################
# ################                           Main Code                    -------------------------------------------
####################################################################################################################



########################################## Basic Seting ------------------------------------------------------------



#set.seed(1234)

#Sample Size
N1<-100
N2<-200
N3<-300
NSize = c(N1, N2, N3) 


##For generating spatial weighting matrix

TrueParaList = list(c(1,1,1),
                    c(5, 10, 20),
                    c(2.5, 2.5, 2.5))

#Range of quantiles
taus = c(0.1, 0.5, 0.9)

######MC times
Krep =1000

#Specify the Spatile weighting matrix
m =3 # network type: 1.DyadW 2.BlockW 3.PowerLawW 4. the simple asymetric matrix
# #if m=4, we also need to specify the number of neighbors
# nn<-2
TruePara = TrueParaList[[m]] 
k=1
N = NSize[k] # Sample size: NSize = c(100, 200, 500)

#specify time period
BN<-10
Time = NSize[1]+BN
#Experiments
exp<-1

#Number of parameters
NP<-13

##IVs
#IVs<-"wx12"
IVs<-"wy23"
#IVs<-"wy2"
#or "wy23";wy2
#Bandwidth
bw<-4 
#We consider 4 bandwidth: HS,3hs;hb,0.6hb, see Koenker and Xiao 2006
##for grid search
step<-0.005
lstep=30

#For generating the individual invariant varianle and time invariant
ZSigma = 0.5^abs(outer(1:5, 1:5, "-"))
FLag = 1


### save results
thetaEst = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2)
Bias = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2) 
#We consider 4 bandwidth: HS,3hs;hb,0.6hb, see Koenker and Xiao 2006
CVprob1 = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2) 
CVprob2 = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2) 
CVprob3 = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2) 
CVprob4 = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2) 

#For storing estimation without IV

thetaEstIVno = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2)
BiasIVno = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2) 
CVprobIVno = rep(list(rep(list(matrix(0, nrow = NP, ncol = Krep)), length(taus))), 2) 


netDensity = rep(0, length(NSize)*3)                                                                             ### store the network density
TNOE = rep(0, length(NSize)*3) 




# ################################## Starting Simulation -----------------------------------------------------


for (t in 1:2)  # distribution t=1:normal; t=2:t(5)
{
  cat("Spatial Weights: ", m, "N: ", N, "Time: ", Time, "Distribution: ", t,"\n")
  
  ### spatial weighting matrix
  if (m==1)
    W = getDyadW(N, TruePara[k])
  if (m==2)
    W = getBlockW(N, TruePara[k])                                                                                      ### Case 2: Block model: block number 5, 10, 20
  if (m==3)
    W = getPowerLawW(N, TruePara[k]) 
  if (m==4)  # the simple asymetric matrix
    W = getAsymW(N,nn)
  
  wma = as.matrix(W)
  
  # distribution: Z normal; T: t distribution
  type = c("Z", "T")[t]
  if (t==1){
    theta = sapply(taus, function(tau) c(qnorm(tau), fphi(tau),fphi1(tau),fphi2(tau),get.gamma.nodal(qnorm(tau)),get.gamma.exog(qnorm(tau))))
   } else{
    theta = sapply(taus, function(tau) c(qt(tau,5), fphi(tau),fphi1(tau),fphi2(tau),get.gamma.nodal(qt(tau,5)),get.gamma.exog(qt(tau,5))))

  }
  
  
  ############# Starting replication
  for (i in 1:Krep)   # simulation replication
  { 
    #Generating of Z(individual invariant) and F (time invariant)
    Z = mvrnorm(n = N, mu = rep(0, nrow(ZSigma)), Sigma = ZSigma)
    Fat = matrix(rnorm((Time+FLag)*2, 0, 1), ncol = 2)
    # quantile levels
    for (q in 1:length(taus))  
    { 
      oldtime<-Sys.time()
      tau = taus[q]
      cat("MC:",i,"Quantile:",tau, "\n")
   
      
      #Generating the dependent variable, the function simu.data
      Ymat = simu.data(W, Time, Z, Fat, type)
      
      dim(Ymat)
      
      #Specifying the ranging of the grid search for the spatial variable
      # There is no general way determining the range for grid search, and multiple tries could be needed for satisfactory results
      #we change the range of grid search to save time
      SR = seq(theta[2, q] - lstep*step, theta[2, q] + lstep*step, step)
      #Estimation, using the function netQuantileIVdir()
      ThetaEstQuant = netQuantileIVdir(Ymat,BN, W, Z, Fat, tau,SR,IVs)  #, SR)
      
      # save results
      thetaEst[[t]][[q]][, i] = ThetaEstQuant$coefficients
      CVprob1[[t]][[q]][, i] = as.numeric(ThetaEstQuant$coefficients - 1.96*ThetaEstQuant$se[,1] <= theta[, q] 
                                          &  theta[, q]  <= ThetaEstQuant$coefficients + 1.96*ThetaEstQuant$se[,1] )
      CVprob2[[t]][[q]][, i] = as.numeric(ThetaEstQuant$coefficients - 1.96*ThetaEstQuant$se[,2] <= theta[, q] 
                                          &  theta[, q]  <= ThetaEstQuant$coefficients + 1.96*ThetaEstQuant$se[,2] )
      CVprob3[[t]][[q]][, i] = as.numeric(ThetaEstQuant$coefficients - 1.96*ThetaEstQuant$se[,3] <= theta[, q] 
                                          &  theta[, q]  <= ThetaEstQuant$coefficients + 1.96*ThetaEstQuant$se[,3] )
      CVprob4[[t]][[q]][, i] = as.numeric(ThetaEstQuant$coefficients - 1.96*ThetaEstQuant$se[,4] <= theta[, q] 
                                          &  theta[, q]  <= ThetaEstQuant$coefficients + 1.96*ThetaEstQuant$se[,4] )
      
      
      #Estimation without using IV
      thetaEstIVno[[t]][[q]][, i] = ThetaEstQuant$estIVno
      CVprobIVno[[t]][[q]][, i] = as.numeric(ThetaEstQuant$coefficients - 1.96*ThetaEstQuant$seIVno <= theta[, q] 
                                             &  theta[, q]  <= ThetaEstQuant$coefficients + 1.96*ThetaEstQuant$seIVno )
      
      
      Bias[[t]][[q]][, i] = ThetaEstQuant$coefficients - theta[, q]
      cat("Bias:",Bias[[t]][[q]][2, i],lstep*step, step,ThetaEstQuant$se[2,1], "\n")
      BiasIVno[[t]][[q]][, i] = ThetaEstQuant$estIVno - theta[, q]
      
    }
    newtime<-Sys.time()
    timeused<-newtime-oldtime
    cat("Time Used:",timeused, "\n")
  }  # i
  
  #Network density
  TNOE[t] = sum(W>0)                    ### the total number of edges
  netDensity[t] = TNOE[t]/(N^2-N)*100   ### the network density
  
}# t

################################## save the estimation result ######################################################

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
RMSE = matrix(0, nrow = 6, ncol = NP)
RMSE[1, ] = sqrt(rowMeans(Bias[[1]][[1]]^2))
RMSE[2, ] = sqrt(rowMeans(Bias[[1]][[2]]^2))
RMSE[3, ] = sqrt(rowMeans(Bias[[1]][[3]]^2))
RMSE[4, ] = sqrt(rowMeans(Bias[[2]][[1]]^2))
RMSE[5, ] = sqrt(rowMeans(Bias[[2]][[2]]^2))
RMSE[6, ] = sqrt(rowMeans(Bias[[2]][[3]]^2))



mCVprob1 = matrix(0, nrow = 6, ncol = NP)
mCVprob2 = matrix(0, nrow = 6, ncol = NP)
mCVprob3 = matrix(0, nrow = 6, ncol = NP)
mCVprob4 = matrix(0, nrow = 6, ncol = NP)

mCVprob1[1, ] = rowMeans(CVprob1[[1]][[1]])
mCVprob1[2, ] = rowMeans(CVprob1[[1]][[2]])
mCVprob1[3, ] = rowMeans(CVprob1[[1]][[3]])
mCVprob1[4, ] = rowMeans(CVprob1[[2]][[1]])
mCVprob1[5, ] = rowMeans(CVprob1[[2]][[2]])
mCVprob1[6, ] = rowMeans(CVprob1[[2]][[3]]) 


mCVprob2[1, ] = rowMeans(CVprob2[[1]][[1]])
mCVprob2[2, ] = rowMeans(CVprob2[[1]][[2]])
mCVprob2[3, ] = rowMeans(CVprob2[[1]][[3]])
mCVprob2[4, ] = rowMeans(CVprob2[[2]][[1]])
mCVprob2[5, ] = rowMeans(CVprob2[[2]][[2]])
mCVprob2[6, ] = rowMeans(CVprob2[[2]][[3]])

mCVprob3[1, ] = rowMeans(CVprob3[[1]][[1]])
mCVprob3[2, ] = rowMeans(CVprob3[[1]][[2]])
mCVprob3[3, ] = rowMeans(CVprob3[[1]][[3]])
mCVprob3[4, ] = rowMeans(CVprob3[[2]][[1]])
mCVprob3[5, ] = rowMeans(CVprob3[[2]][[2]])
mCVprob3[6, ] = rowMeans(CVprob3[[2]][[3]])

mCVprob4[1, ] = rowMeans(CVprob4[[1]][[1]])
mCVprob4[2, ] = rowMeans(CVprob4[[1]][[2]])
mCVprob4[3, ] = rowMeans(CVprob4[[1]][[3]])
mCVprob4[4, ] = rowMeans(CVprob4[[2]][[1]])
mCVprob4[5, ] = rowMeans(CVprob4[[2]][[2]])
mCVprob4[6, ] = rowMeans(CVprob4[[2]][[3]])



mBias = matrix(0, nrow = 6, ncol = NP)
mBias[1, ] = rowMeans(Bias[[1]][[1]])
mBias[2, ] = rowMeans(Bias[[1]][[2]])
mBias[3, ] = rowMeans(Bias[[1]][[3]])
mBias[4, ] = rowMeans(Bias[[2]][[1]])
mBias[5, ] = rowMeans(Bias[[2]][[2]])
mBias[6, ] = rowMeans(Bias[[2]][[3]]) 


#Without using IV
RMSEIVno = matrix(0, nrow = 6, ncol = NP)
RMSEIVno[1, ] = sqrt(rowMeans(BiasIVno[[1]][[1]]^2))
RMSEIVno[2, ] = sqrt(rowMeans(BiasIVno[[1]][[2]]^2))
RMSEIVno[3, ] = sqrt(rowMeans(BiasIVno[[1]][[3]]^2))
RMSEIVno[4, ] = sqrt(rowMeans(BiasIVno[[2]][[1]]^2))
RMSEIVno[5, ] = sqrt(rowMeans(BiasIVno[[2]][[2]]^2))
RMSEIVno[6, ] = sqrt(rowMeans(BiasIVno[[2]][[3]]^2))


mCVprobIVno = matrix(0, nrow = 6, ncol = NP)

mCVprobIVno[1, ] = rowMeans(CVprobIVno[[1]][[1]])
mCVprobIVno[2, ] = rowMeans(CVprobIVno[[1]][[2]])
mCVprobIVno[3, ] = rowMeans(CVprobIVno[[1]][[3]])
mCVprobIVno[4, ] = rowMeans(CVprobIVno[[2]][[1]])
mCVprobIVno[5, ] = rowMeans(CVprobIVno[[2]][[2]])
mCVprobIVno[6, ] = rowMeans(CVprobIVno[[2]][[3]]) 




mBiasIVno = matrix(0, nrow = 6, ncol = NP)
mBiasIVno[1, ] = rowMeans(BiasIVno[[1]][[1]])
mBiasIVno[2, ] = rowMeans(BiasIVno[[1]][[2]])
mBiasIVno[3, ] = rowMeans(BiasIVno[[1]][[3]])
mBiasIVno[4, ] = rowMeans(BiasIVno[[2]][[1]])
mBiasIVno[5, ] = rowMeans(BiasIVno[[2]][[2]])
mBiasIVno[6, ] = rowMeans(BiasIVno[[2]][[3]]) 




############################################### Make Latex Table --------------------------------------------------------


##Make Latex table for RMSE

latex = rbind(paste("\\hline & & \\multicolumn{13}{c}{$ N = ", N, ", Network Density = ",  specify_decimal(netDensity[1], 2), "%; ",
                    specify_decimal(netDensity[2], 2), "% $}"),
              paste("\\hline & & \\multicolumn{13}{c}{RMSE}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*RMSE[1, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*RMSE[2, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*RMSE[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*RMSE[4, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*RMSE[5, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*RMSE[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Converge Probability (1)}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*mCVprob1[1, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob1[2, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob1[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*mCVprob1[4, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob1[5, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob1[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Converge Probability (2)}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*mCVprob2[1, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob2[2, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob2[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*mCVprob2[4, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob2[5, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob2[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Converge Probability (3)}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*mCVprob3[1, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob3[2, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob3[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*mCVprob3[4, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob3[5, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob3[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Converge Probability (4)}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*mCVprob4[1, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob4[2, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob4[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*mCVprob4[4, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprob4[5, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprob4[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Bias}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*mBias[1, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*mBias[2, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*mBias[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*mBias[4, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*mBias[5, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*mBias[6, ], 2), collapse = " & ")), 
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Without using IV}"),
              paste("\\hline & & \\multicolumn{13}{c}{RMSE}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*RMSEIVno[1, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*RMSEIVno[2, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*RMSEIVno[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*RMSEIVno[4, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*RMSEIVno[5, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*RMSEIVno[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Converge Probability }"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*mCVprobIVno[1, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprobIVno[2, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprobIVno[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*mCVprobIVno[4, ], 2), collapse = " & ")),
              paste(" &  0.5 &", paste(specify_decimal(100*mCVprobIVno[5, ], 2), collapse = " & ")),
              paste(" &  0.9 &", paste(specify_decimal(100*mCVprobIVno[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"),
              paste("\\hline & & \\multicolumn{13}{c}{Bias}"),
              paste("$N(0,1)$ &  0.1 &", paste(specify_decimal(100*mBiasIVno[1, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*mBiasIVno[2, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*mBiasIVno[3, ], 2), collapse = " & ")),
              paste("$t(5)$ &  0.1 &", paste(specify_decimal(100*mBiasIVno[4, ], 2), collapse = " & ")),
              paste("    &  0.5 &", paste(specify_decimal(100*mBiasIVno[5, ], 2), collapse = " & ")),
              paste("    &  0.9 &", paste(specify_decimal(100*mBiasIVno[6, ], 2), collapse = " & ")),
              paste("\\cline{2-16}"))

#Save the latex table
print(xtable(latex, type = "latex"), file= sprintf("results/Bias, RMSE and CP, Weights=%d (Wpara=%0.1f),IV=%s,N=%d,T=%d,MC=%d,bandwidth=%d,SR=(%0.4f,%d),IVdir(%d).tex",m,TruePara[k],IVs,N,Time,Krep,bw,step,lstep,exp))


cat("Spatial Weights: ", m, "N: ", N, "Time: ", Time, "Distribution: ", t,"IV",IVs, "netDensity:", netDensity,"\n")
# mBiaslatex
# Convergencelatex    
# RMSElatex

result_all<-rbind(mBias,RMSE,mCVprob1,mCVprob2,mCVprob3,mCVprob4,mBiasIVno,RMSEIVno,mCVprobIVno)
write.csv(result_all,file= sprintf("results/Bias, RMSE and CP, Weights=%d (Wpara=%0.1f),IV=%s,N=%d,T=%d,MC=%d,bandwidth=%d,SR=(%0.4f,%d),IVdir(%d).csv",m,TruePara[k],IVs,N,Time,Krep,bw,step,lstep,exp))




