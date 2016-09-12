vcov.cglm <- function(object, type = "model-based") {  
    model <- object
    if (!inherits(model, "cglm")) stop("model must be an object of class \"cglm\"")
    if (!(type %in% c("model-based", "robust"))) {
        stop("type must be either \"model-based\" or \"robust\"")
    }
    X <- model$x
    y <- model$y
    totals <- rowSums(y)
    n <- nrow(y)
    d <- ncol(y)
    coefs <- coef(model)
    np <- nrow(coefs)
    fullContrasts <- as.matrix(model$coef.contrasts)
    ##
    ##  Model-based vcov matrix uses the model-based estimate of "Sigma":
    ##
    if (type == "model-based"){
        Sig <- model$Sigma
        result <- kronecker(crossprod(fullContrasts, Sig) %*% fullContrasts, solve(crossprod(X)))
    }
    ##
    if (type == "robust"){
        ##
        if (model$method %in% c("logy.fit", "multinom.fit")) { 
            stop("robust standard errors are not implemented for this model")
        }
        ##  The bread in the Sandwich estimator is model-based vcov with Sigma the identity matrix
bread <- kronecker(crossprod(fullContrasts, fullContrasts), solve(crossprod(X)))
        ##  Next compute the "ham": 
        p <- fitted(model)
        p <- p / totals
        y <- y / totals
        res <- (y - p)
        S <-                              ## array of row-residual SSP matrices
            lapply ((1:n),  function(i) {
                res.i <- res[i,]
                tcrossprod(res.i)
            })
        S <- bdiag(S)
        PP <-                              ## "multinomial" variance-covariance matrices
            lapply ((1:n),  function(i) {
                p.i <- p[i,]
                #        crossprod(theContrasts, diag(p.i) - outer(p.i, p.i)) 
                diag(p.i) - tcrossprod(p.i)
            })
        Xexp <- kronecker(X, diag(d))      ## expanded to all D components
        PPexp <- bdiag(PP)
        M <- diag(d) - 1/d
        Vinv <- lapply(1:n, function(i) {
            p.i <- p[i,]
            Pinv <- diag(1/p.i)
            M %*% Pinv %*% M %*% Pinv %*% M
        })
        Vinv <- bdiag(Vinv)
        D <- PPexp %*% Xexp
        VinvD <- Vinv %*% D
        ham <- as.matrix(Matrix::crossprod(VinvD, S) %*% VinvD)
        perm <- as.vector(outer(d * (0:(np-1)), 1:d, FUN = "+"))  
        ham <- ham[perm, perm]  ## corrects the ordering of rows and columns
        ##
        result <- bread %*% ham %*% bread    ## the sandwich formula
        ##  Now extract the vcov relevant to the specified set of contrasts:
        cc <- kronecker(fullContrasts, diag(np))
        result <- crossprod(cc, result) %*% cc
    }
    ##
    ##  Finally, tidy up the row and column names in the result:
    coefnames <- dimnames(coefs)
    coefnames <- paste(rep(coefnames[[1]], ncol(coefs)), 
                       rep(coefnames[[2]], rep(nrow(coefs), ncol(coefs))), 
                       sep = "_")
    rownames(result) <- colnames(result) <- coefnames
 
   return(result)
} 
