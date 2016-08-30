gw.fit <- function(X, Y, weights = rep(1, nobs), offset = rep(0, nobs), 
                     maxit = 20, epsilon = 1e-8, trace = FALSE, 
                     intercept = TRUE, br = FALSE, checkAliasing = FALSE, qr = FALSE){
  
  if (br) warning("br=TRUE has been ignored (it is not implemented for this method)")
  nobs <- nrow(Y)
  totals <- rowSums(Y)
  log.totals <- log(totals)
  D <- ncol(Y)  
  
  ###  Next function does one iteration of the hybrid quasi-likelihood algorithm
  hybrid.update <- function(fit){
    eta <- fitted(fit)
    mu <- exp(eta)
    fitted.totals <- rowSums(mu)
    mu <- mu * totals / fitted.totals
    z <- Y/mu - 1
    adj.dep.var <- eta + z - rowMeans(z)
    fit <- lm(adj.dep.var ~ X - 1, weights = weights)
    return(fit)
  }
 
  ###  Initialize the iterative algorithm
  Yadj <- Y + 0.01  ## to avoid log(0); only used just here!  (not later)
  logY <- log(Yadj) ## we'll use this as working dependent variate; just here
  fit <- lm(logY ~ -1 + X, weights = weights)
  iteration <- 0
  converged <- FALSE
  
  ###  Do the hybrid iterations
  while((iteration < maxit) && !converged){
    iteration <- iteration + 1
    oldfit <- fit
    fit <- hybrid.update(fit)
    if (trace) {
      cat("Iteration ", iteration, "\n")
      print(fit)
    }
    converged <- all((abs(fitted(fit) - fitted(oldfit))) < epsilon)
  }
  
  ###  Tidy up and finish cglm.fit 
  rownames(fit$coefficients) <- colnames(X)
  fit$fitted.values <- exp(fit$fitted.values)
  fit$fitted.values <- totals * fit$fitted.values / rowSums(fit$fitted.values)
  fit$residuals <- Y / fit$fitted.values - 1    
  fit$residuals <- fit$residuals - rowMeans(fit$residuals)  ## arbitrary, but symmetric in the components
  ##
  fit$Sigma <- crossprod(fit$residuals) / fit$df.residual   ## estimates (C_*)' Sigma (C_*)
  ##
  fit$iter <- iteration
  fit$converged <- converged
  if (!converged) warning("not converged")
  return(fit)
}
