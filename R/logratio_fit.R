logratio_fit <- function(X, Y, weights = rep(1, nobs), offset = rep(0, nobs),
                     maxit = NULL, epsilon = NULL, trace = FALSE,
                     intercept = TRUE){
  nobs <- nrow(Y)
  if (all(Y > 0)) logY <- log(Y)
  else stop("logy.fit can only be used when all y-values are positive")
  fit <- lm(logY ~ X - 1, weights = weights)
  fit$linear.predictors <- fit$fitted.values
  fitted <- exp(fit$fitted.values)
  ##
  fit$fitted.values <- fitted * rowSums(Y) / rowSums(fitted)  ## ensures fitted row sums agree with those of Y
  ##
  fit$residuals <- fit$residuals - rowMeans(fit$residuals)    ## arbitrary, but symmetric in the components
  ##
  null.fit <- colSums(weights * logY) / sum(weights)
  null.res <- t(logY) - null.fit
  null.res <- null.res - colMeans(null.res)
  fit$SS.null <- sum(weights * (null.res ^ 2))
  ##
  fit$Sigma <- crossprod(fit$residuals) / fit$df.residual     ## estimates (C_*)' Sigma (C_*)
  ##
  rownames(fit$coefficients) <- colnames(X)
  colnames(fit$coefficients) <- colnames(Y)
  fit$iter <- 1
  fit$converged <- TRUE
  return(fit)
}
