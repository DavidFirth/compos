cglm <- function (formula, data, weights, subset,
                    na.action, start = NULL,
                    maxit = 100, epsilon = 1e-8, trace = FALSE,
                    model = TRUE, method = "gw.fit",
                    y = TRUE, contrasts = NULL, ref = "mean",
                  y.totals = NULL, br = FALSE, checkAliasing = FALSE, qr = FALSE, ...)
    ## IK NOTE 18/06/2016: The chekAliasing and qr areguments have
    ## been temporary added here. They work with multinom.fit and are
    ## there to control whether a qr decompositions needs to be
    ## computed on exit (useful for vcov methods) and/or whether
    ## aliasing checks need to made. The addition of these two
    ## arguments may affect the other methods, if they are specific in
    ## the arguments they take.
    ##
    ##  The y.totals arg is there to allow NAs to be treated properly
    ##  (in due course, once I have programmed that)
{
    call <- match.call()
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame"))
        return(mf)
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    case.names <- rownames(Y)
    X <- if (!is.empty.model(mt)) {
        model.matrix(mt, mf, contrasts)
    } else matrix(, NROW(Y), 0L)
    rownames(X) <- case.names

    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")

###  Set up the Y and X matrices
    if (!is.matrix(Y)) stop("The y-variable must be a matrix")
    DD <- ncol(Y)
    n <- nrow(Y)
    if (is.null(colnames(Y))) {
      colnames(Y) <- paste("c", 1:DD, sep = "")
    }
    response.names <- colnames(Y)
    ones <- rep(1, DD); names(ones) <- response.names
    rsY <- rowSums(Y)
##  Do a check on the specified y.totals here, and report the largest errors
##  -- or stop with an error if y.totals is not positive or a vector of length 1 or n
##  The specific error checking that's done next should be avoided, though
      if (is.numeric(y.totals) && (y.totals > .Machine$double.eps ^ 0.5)) {
	    if (!isTRUE(all.equal(rsY / y.totals, rep(1, n), check.attributes = FALSE))) {
	        stop("The specified y.totals value seems to be wrong")
	    }
      Y <- Y / y.totals
    }
#    else if (is.null(y.totals)) {
#	    Y <- Y / rsY
#    } else stop("y.totals must be either a positive number or NULL")

offset <- rep(0, n)
##  Still need to handle offset in the below.  A vector of length nrow(X).

    fit <- eval(call(if (is.function(method)) "method" else method,
                  X = X, Y = Y, weights = weights, offset = offset,
                  maxit = maxit, epsilon = epsilon, trace = trace >
                  0L, br = br, checkAliasing = checkAliasing, qr = qr))
    coefs <- coef(fit)
    fit$coef.contrasts <- cmatrix <- cglmContrastMatrix(coefs, ref)
    fit$coefficients <- as.matrix(coefs %*% cmatrix)
    fit$ref <- ref
    fit$call <- call
    fit$y <- Y
    fit$x <- X
    fit$prior.weights <- if (is.null(weights)) rep(1, n) else weights
    names(fit$prior.weights) <- case.names
    fit$linear.predictors <- log(fit$fitted.values)
    fit$na.action = attr(mf, "na.action")
    fit$method <- method
    fit$deviance <- deviance.cglm(fit)
    if (model) fit$model <- mf
    class(fit) <- c("cglm", "mlm", "lm")
    return(fit)
}



