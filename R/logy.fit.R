logy.fit <- function(X, Y, weights = rep(1, nobs), offset = rep(0, nobs), 
                     maxit = NULL, epsilon = NULL, trace = FALSE, 
                     intercept = TRUE, br = FALSE, checkAliasing = FALSE, qr = FALSE){
  if (br) warning("br=TRUE has been ignored (it is not implemented for this method)")
  nobs <- nrow(Y)
  if (all(Y > 0)) logY <- log(Y) 
  else stop("logy.fit can only be used when all y-values are positive")
  fit <- lm(logY ~ X - 1, weights = weights)
  fitted <- exp(fit$fitted.values)
  ##
  fit$fitted.values <- fitted * rowSums(Y) / rowSums(fitted)  ## ensures fitted row sums agree with those of Y
  ##
  fit$residuals <- fit$residuals - rowMeans(fit$residuals)    ## arbitrary, but symmetric in the components
  ##
  fit$Sigma <- crossprod(fit$residuals) / fit$df.residual     ## estimates (C_*)' Sigma (C_*)
  ##
  rownames(fit$coefficients) <- colnames(X)
  colnames(fit$coefficients) <- colnames(Y) 
  fit$iter <- 1
  fit$converged <- TRUE
  return(fit)
}
