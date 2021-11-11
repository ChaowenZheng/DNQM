
# network estimation
netQuantile<-function(Ymat, W, Z, M, tau = 0.5, Fat = NULL, A = NULL)
{
  WYmat = W%*%Ymat
  Time = ncol(Ymat) - 1
  
  # endogenous variable D
  D = WYmat[, -1]
  D = as.vector(D)
  
  if (is.null(Fat))
  {
    # exogenous variable X
    X = cbind(as.vector(WYmat[, -ncol(Ymat)]),       # gamma2
              as.vector( Ymat[, -ncol(Ymat)]),       # gamma3  
              do.call("rbind", rep(list(Z), Time)))  # alpha1-alpha4  4 "SIZE", "BM", "Cash", "PE" 

  }
  else
  {
    f = cbind(Fat[2:nrow(Fat), 1], Fat[1:nrow(Fat)-1, 1],
              Fat[2:nrow(Fat), 2], Fat[1:nrow(Fat)-1, 2])
    
    Fatlist = lapply(1:4, function(i) { Ft = t(f[, i]%x%matrix(1, 1, nrow(W)))
    return(as.vector(Ft)) })
    
    # exogenous variable X
    X = cbind(as.vector(WYmat[, -ncol(Ymat)]),   # gamma2
              as.vector( Ymat[, -ncol(Ymat)]),    # gamma3  
              do.call("rbind", rep(list(Z), Time)), # alpha1-alpha4 4 "SIZE", "BM", "Cash", "PE"
              do.call("cbind", Fatlist))  # beta1-beta4 4 "VIX", "Rm - Rf", "SMB", "HML"
  }
  
  # M is IV, XIV including X and IV
  XIV = cbind(X, M)
  
  # dependent variable
  Yvec = as.vector(Ymat[, -1])
  
  ## Step 1 obtain the range of A or alpha parameter
  if (is.null(A))
  {
    ## the initial value 
    #  independent variable X without IV
    if (is.null(Fat))
    {
      # exogenous variable X
      XIVnon = cbind(as.vector(WYmat[, -1]),               # gamma1
                     as.vector(WYmat[, -ncol(Ymat)]),      # gamma2 lagged network effects
                     as.vector( Ymat[, -ncol(Ymat)]),      # gamma3  
                     do.call("rbind", rep(list(Z), Time))) # alpha1-alpha4  4 "SIZE", "BM", "Cash", "PE"  
      
    }
    else
    {
      # exogenous variable X
      XIVnon = cbind(as.vector(WYmat[, -1]),               # gamma1
                     as.vector(WYmat[, -ncol(Ymat)]),      # gamma2
                     as.vector( Ymat[, -ncol(Ymat)]),      # gamma3  
                     do.call("rbind", rep(list(Z), Time)), # alpha1-alpha4 4  "SIZE", "BM", "Cash", "PE"
                     do.call("cbind", Fatlist))            # beta1-beta4 4    "VIX", "Rm - Rf", "SMB", "HML"
    }
    
    qrnonIV = rq(Yvec ~ XIVnon, tau)
    qrnonIV
    qrnonIV_coef = summary.rq(qrnonIV)$coefficients 
    gamma1_IVnon = qrnonIV_coef[2, 1]
    cat("\n ", gamma1_IVnon)
    std_IVnon = qrnonIV_coef[2, 2]
    alpha = seq(gamma1_IVnon - 15*std_IVnon, 
                gamma1_IVnon + 15*std_IVnon, 3*std_IVnon)
  }
  else
  {
    alpha = A
  }
  

  ## Step 2 get alpha by the minimum value of IV
  gnorm = lapply(alpha, function(t)
  {         # t = alpha[1]
    cat(" ", t) 
    dat2   = data.frame(Y = Yvec - t*D, XIV)
    resrq2 = rq(Y~., tau, data = dat2)
    # infer_coef = summary.rq(resrq, se = "boot", hs = F)$coefficients   # boot
    # infer_coef
    # g = norm(as.matrix(tail(resrq2$coefficients, ncol(M))), "F")           # endogenous D dimension 1
    g = norm(as.matrix(tail(resrq2$coefficients, ncol(M))*M), "F")
    # g = norm(as.matrix(tail(infer_coef[, 1], 1)), "F")           # endogenous D dimension 1
    return(g)
  })
  gnorm  = do.call("rbind", gnorm)
  gamma1 = alpha[which(gnorm == min(gnorm))]
  cat("\n ", gamma1)
  
  ## Step 3 plug in alpha and get other estimators
  dat   = data.frame(Y = Yvec - gamma1*D, X)
  resrq = rq(Y~., tau, data = dat)
  # infer_coef = summary.rq(resrq, se = "boot", hs = F)$coefficients
  infer_coef = resrq$coefficients

  # combine the coefficients with gamma1
  esti_coef = c(infer_coef[1], gamma1, infer_coef[2:(ncol(X)+1)])
  cat("\n ", esti_coef)
  
  
  ### 
  dat3   = data.frame(Y = Yvec - gamma1*D, XIV)
  resrq3 = rq(Y~., tau, data = dat3)
  R = as.matrix(M) %*% as.matrix(tail(resrq3$coefficients, ncol(M)))
  # R = M[, 1]
  
  # calculate the se, t value using the coefficents
  bout = vciqr(bhat = esti_coef, y = Yvec, d = D, x =  X, z = R, tau = tau)

  # output 
  return(bout)
}


## compare without simutanous term
lagnetQuantile<-function(Ymat, W, Z, tau = 0.5, Fat = NULL)                                                                                  ### OLS estimation for theta: eq (2.8)
{
  WYmat = W%*%Ymat                                                                                             ### obtain WY
  Time = ncol(Ymat)-1
  
  if (is.null(Fat))
  {
    # exogenous variable X
    X = cbind(as.vector(WYmat[, -ncol(Ymat)]),       # gamma2
              as.vector( Ymat[, -ncol(Ymat)]),       # gamma3  
              do.call("rbind", rep(list(Z), Time))) # alpha1-alpha4  4 "SIZE", "BM", "Cash", "PE" 
    
  }
  else
  {
    f = cbind(Fat[2:nrow(Fat), 1], Fat[1:nrow(Fat)-1, 1],
              Fat[2:nrow(Fat), 2], Fat[1:nrow(Fat)-1, 2])
    
    Fatlist = lapply(1:4, function(i) { Ft = t(f[, i] %x% matrix(1, 1, nrow(W)))
    return(as.vector(Ft)) })
    
    # exogenous variable X
    X = cbind(as.vector(WYmat[, -ncol(Ymat)]),    # gamma2
              as.vector( Ymat[, -ncol(Ymat)]),    # gamma3  
              do.call("rbind", rep(list(Z), Time)), # alpha1-alpha4 4 "SIZE", "BM", "Cash", "PE"
              do.call("cbind", Fatlist))  # beta1-beta4 4 "VIX", "Rm - Rf", "SMB", "HML"
  }
  
  Yvec = as.vector(Ymat[, -1])                                                                                  ### the response vector
  dat = data.frame(Y = Yvec, X)
  resrq = rq(Y~., tau , data = dat)
  infer_coef = summary.rq(resrq, se = "nid", hs = F)$coefficients
  lagresid = as.matrix(resid(resrq))
  
  return(list(coefficients = infer_coef, res = lagresid, tau = tau, model = resrq))
}



# b = vciqr(bhat = esti_coef, y = Yvec, d = D, x =  X, z = Dhat, tau = 0.9)
# # Compute standard errors for IQR when the equation is just identified
# # using the method proposed in Powell (1986).  This implementation uses
# % a Gaussian kernel and a simple rule of thumb bandwidth.
# % [b,vc,J] = vciqr(bhat,y,x,tau)
# % Inputs -
#   %   bhat - estimated coefficients from QR
# %   y - dependent variable
# %   d - RHS endogenous variable
# %   x - covariates (If there are no covariates, pass x = [].)
# %   z - instruments
# %   tau - quantile index
# % Outputs -
#   %   b - estimated coefficients with standard errors
# %   vc - covariance matrix of b
# %   J - matrix in asymptotic variance formula (J'SJ) used in process testing
vciqr <- function(bhat, y, d, x, z, tau = 0.9)
{
  
n = length(y)	    # Number of observations
# x = cbind(1, x)  # Add constant term
xd = cbind(1, d, x)
k = ncol(xd)       #   Number of regressors

vc = matrix(0, nrow = k, ncol = k)
b = matrix(0, nrow = k, ncol = 4)
xz = cbind(1, z, x)
S = (1/n)*t(xz)%*%xz
                                               
e = y - xd%*%bhat     # Generate residuals

h = 1.364*((2*sqrt(pi))^(-1/5))*sqrt(var(e))*(n^(-1/5))  # Calculate bandwidth using Silverman's rule of thumb
J = as.vector(1/(n*h))*(t((dnorm(as.vector(1/h)*e)%*%rep(1, ncol(xd)))*xd)%*%xz)

vc = (1/n)*(tau-tau^2)*solve(t(J))%*%S%*%solve(J)

b[, 1] = bhat
b[, 2]  = (sqrt(diag(vc)))
b[, 3] = bhat / (sqrt(diag(vc)))
b[, 4] = 2*pnorm(-abs(bhat / (sqrt(diag(vc)))))

cnames = c( "Value", "Std. Error", "t value", "Pr(>|t|)")
if (k < 9){
  rnames = c("intercept", "gamma1", "gamma2", "gamma3", "SIZE", "BM", "Cash", "PE")
  }
else{
  rnames = c("intercept", "gamma1", "gamma2", "gamma3", "SIZE", "BM", "Cash", "PE", "VIX", "Rm - Rf", "SMB", "HML")
} 
dimnames(b) = list(rnames, cnames)
                           
J = solve(J)

return(list(b = b, vc = vc, J = J, res = e, tau = tau))
}


goodfit <- function(resid, resid_nl, tau){
  # minimum sum of deviations
  V1 <- resid * (tau - (resid < 0))
  V1 <- sum(V1, na.rm = T) 
  
  # null sum of deviations
  V0 <- resid_nl * (tau - (resid_nl < 0))
  V0 <- sum(V0, na.rm = T) 
  
  # explained deviance
  out <- 1 - V1/V0
  
  # # exceptions for output
  # if(any(c(Inf, -Inf) %in% out)) out <- NA
  # if(V1 > V0) out <- NA
  
  return(out)
}





Plotnetivrq <- function (x, parm = NULL, level = 0.9, mfrow = NULL, 
          mar = NULL, ylim = NULL, main = NULL, col = gray(c(0, 0.75)), 
          border = NULL, lcol = 2, lty = 1:2, cex = 0.5, pch = 20, 
          type = "b", xlab = "", ylab = "") 
{
  zalpha <- qnorm(1 - (1 - level)/2)
  taus <- sapply(x, function(x) x$tau)
  
  cf = lapply(x, function(x) x$b)
  if (ncol(cf[[1]]) == 4) {
    for (i in 1:length(cf)) {
      cfi <- cf[[i]]
      cfi <- cbind(cfi[, 1], cfi[, 1] - cfi[, 2] * zalpha, 
                   cfi[, 1] + cfi[, 2] * zalpha)
      colnames(cfi) <- c("coefficients", "lower bd", "upper bd")
      cf[[i]] <- cfi
    }
  }
  if (ncol(cf[[1]]) != 3) 
    stop("summary.rqs components have wrong dimension")
  if (is.null(parm)) 
    parm <- rownames(cf[[1]])
  if (is.numeric(parm)) 
    parm <- rownames(cf[[1]])[parm]

  cf <- lapply(cf, function(x) x[parm, , drop = FALSE])
  names(cf) <- paste("tau=", taus)
  
  mfrow_orig <- par("mfrow")
  mar_orig <- par("mar")
  if (is.null(mfrow)) 
    mfrow <- n2mfrow(length(parm))
  if (is.null(mar)) 
    mar <- c(3.1, 3.1, 3.1, 1.6)
  par(mfrow = c(4,3), mar = mar)
  col <- rep(col, length.out = 2)
  lty <- rep(lty, length.out = 2)
  if (is.null(border)) 
    border <- col[2]
  if (is.null(main)) 
    main <- parm

  main <- rep(main, length.out = length(parm))
  xlab <- rep(xlab, length.out = length(parm))
  ylab <- rep(ylab, length.out = length(parm))
  ylim0 <- ylim
  for (i in seq(along = parm)) {
    b <- t(sapply(seq(along = cf), function(tau) cf[[tau]][i, ]))
    if (is.null(ylim))  {
      ylim <- range(b[, 2], b[, 3])
    }
    plot(rep(taus, 2), c(b[, 2], b[, 3]), type = "n", ylim = ylim, 
         xlab = xlab[i], ylab = ylab[i], main = main[i])
    polygon(c(taus, rev(taus)), c(b[, 2], rev(b[, 3])), col = col[2], 
            border = border)
    points(taus, b[, 1], cex = cex, pch = pch, type = type, 
           col = col[1])
    abline(h = 0, col = gray(0.3))
    ylim <- ylim0
  }
  par(mfrow = mfrow_orig, mar = mar_orig)
  x <- cf
  invisible(structure(as.vector(unlist(x)), .Dim = c(dim(x[[1]]), length(x)), .Dimnames = list(rownames(x[[1]]), colnames(x[[1]]), 
                                                                                               names(x))))
}

