cglmContrastMatrix <- function(coef.matrix, ref = NULL) {
  require(Matrix)
  D <- ncol(coef.matrix)
  cnames <- colnames(coef.matrix)
  if (ref == "mean") {  ## centered contrasts are made, in this case
    cmatrix <- Diagonal(D) - Matrix(1, D, D) / D
    rownames(cmatrix) <- colnames(cmatrix) <- cnames
    return(cmatrix)
  }
  else {    ## if not NULL, ref must be a valid column index for coef.matrix
    if (!(ref %in% 1:D) && !(ref %in% cnames)) {
      stop("ref must either be a valid component of the response or else \"mean\"")
    } 
    cmatrix <- Matrix(0, D, D, sparse = TRUE)
    rownames(cmatrix) <- colnames(cmatrix) <- cnames
    cmatrix[ref, -ref] <- -1
    cmatrix[-ref, -ref] <- Diagonal(D - 1) + cmatrix[-ref, -ref]
    return(cmatrix)
  }
}
