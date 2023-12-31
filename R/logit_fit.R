logit_fit <- function(X, Y, weights = rep(1, nobs), offset = rep(0, nobs),
                     maxit = 20, epsilon = 1e-8, trace = FALSE,
                     intercept = TRUE, checkAliasing = FALSE, qr = FALSE){

  nobs <- nrow(Y)
  totals <- rowSums(Y)
  log.totals <- log(totals)
  D <- ncol(Y)

  ###  Next function does one iteration of the iterative least-squares algorithm
  update <- function(fit){
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
    fit <- update(fit)
    if (trace) {
      cat("Iteration ", iteration, "\n")
      print(fit)
    }
    converged <- all((abs(fitted(fit) - fitted(oldfit))) < epsilon)
  }

  ###  Tidy up and finish gw.fit
  rownames(fit$coefficients) <- colnames(X)
  fit$linear.predictors <- fit$fitted.values
  fit$fitted.values <- exp(fit$fitted.values)
  fit$fitted.values <- totals * fit$fitted.values / rowSums(fit$fitted.values)
  fit$residuals <- Y / fit$fitted.values - 1
  fit$residuals <- fit$residuals - rowMeans(fit$residuals)  ## arbitrary, but symmetric in the components
  ##
  null.fit <- colSums(weights * Y) / sum(weights)
  null.res <- t(Y) / null.fit - 1
  null.res <- null.res - colMeans(null.res)
  fit$SS.null <- sum(weights * (null.res ^ 2))
  ##
  fit$Sigma <- crossprod(fit$residuals) / fit$df.residual   ## estimates (C_*)' Sigma (C_*)
  ##
  fit$iter <- iteration
  fit$converged <- converged
  if (!converged) warning("not converged")
  return(fit)
}
