summary.cglm <- function (object, digits = max(3, getOption("digits")-3), vcov_type = "model-based", ...) 
{
    call <- object$call
    ref <- object$ref
    coef <- coef(object)
    refname <- colnames(coef)[ref]
    p <- nrow(coef)
    D <- ncol(coef)
    resid <- residuals(object)
    sterrs <- matrix(sqrt(diag(vcov(object, type = vcov_type))), p, D)[, -ref, drop = FALSE]
    coef <- coef[, -ref, drop = FALSE]
    dimnames(sterrs) <- dimnames(coef)
    ynames <- colnames(coef)
    coeflist <- setNames(vector("list", D-1), paste("  ", ynames, "/", refname))
    for (i in seq(ynames)) {
        coef.se <- cbind(coef[, i], sterrs[, i])
        colnames(coef.se) <- c("Estimate", "St. err")
        rownames(coef.se) <- rownames(coef)
        coeflist[[i]] <- signif(coef.se, digits = digits)
    }
    class(coeflist) <- "listof"
    if (object$df.residual > 5L) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (length(dim(resid)) == 2L) 
            structure(apply(t(resid), 1L, quantile), 
                      dimnames = list(nam, dimnames(resid)[[2L]]))
        else {
            zz <- zapsmall(quantile(resid), digits + 1L)
            structure(zz, names = nam)
        }
    }
    value <- list(
        Call = call,
        Residuals = signif(rq, digits),
        Coefficients = coeflist
    )
    class(value) <- "listof"
    value
}
