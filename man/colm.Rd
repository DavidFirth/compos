\name{colm}
\alias{colm}
\title{
Specify and Fit Multivariate Regression Models to Compositions
}
\description{
This function works in much the same way as \code{stats::lm} or \code{stats::glm}.  The
models currently available are Aitchison-type multivariate linear
regression based on logratio transformation of the responses, and the
generalized Wedderburn logit model as developed in Firth and Sammut (2023).
}
\usage{
colm(formula, data, weights, subset, na.action = na.exclude, start = NULL, maxit = 100, 
   epsilon = 1e-08, trace = FALSE, method = "logit_fit", contrasts = NULL, ref = 1, 
   y.totals = NULL, qr = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{
  a model formula.  See \code{\link{glm}} for more details.
}
\item{data}{
  an optional data frame, list or environment as for \code{\link{glm}}.  If not found in \code{data}, the variables for the model are taken from environment(formula), typically the environment from which colm is called.
}
\item{weights}{
  an optional vector of 'prior weights' to be used in the fitting process. Should be NULL or a numeric vector.
}
\item{subset}{
  an optional vector specifying a subset of observations to be used in the fitting process.
}
\item{na.action}{
  must be either \code{\link{na.omit}} or \code{\link{na.exclude}}.
}
\item{start}{
  a matrix of starting values for the parameters in the linear predictor.
}
\item{maxit}{
  the maximum number of least squares iterations to perform.
}
\item{epsilon}{
  the convergence tolerance.  ?? say how it is defined
}
\item{trace}{
  logical, whether progress of iterations should be reported.
}
\item{contrasts}{
  an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix}}.
}
\item{ref}{
  integer, the response reference component to use.
}
\item{y.totals}{
  numeric, either a scalar or a vector whose length is the sample size.  Not currently used.
}
\item{method}{
  either \code{"logit_fit"} (for generalized Wedderburn logit model) or
  \code{"logratio_fit"} (for Aitchison-type multivariate linear model
  after logratio transformation).
}
\item{qr}{
  logical.  If ‘TRUE’ the corresponding component of the fit (the QR matrix) is returned.
}
\item{\dots}{
  arguments to pass to other functions.
}
}
\value{
A S3 model object of class \code{"colm"}, as a list with the following components:
  \item{coefficients}{a matrix, the parameter estimates in the fitted model.}
  \item{residuals}{the matrix of model residuals.}
  \item{effects}{??}
  \item{rank}{the numeric rank of the fitted model.}
  \item{fitted.values}{fitted values from the model.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
  \item{call}{the matched call.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{model}{the model frame.}
  \item{Sigma}{the estimate of the residual coefficient of variation and covariation matrix.}
  \item{iter}{number of multivariate least squares iterations performed.}
  \item{converged}{logical, whether or not the iterative algorithm has met the specified convergence criterion.}
  \item{coef.contrasts}{the matrix of contrasts implied by the value of \code{ref}}
  \item{ref}{the reference component used.}
  \item{y}{matrix whose rows are the compositional responses for the model.}
  \item{x}{the model matrix (design matrix).}
  \item{prior.weights}{the weights initially supplied, or a vector of 1s if none were.}
  \item{linear.predictors}{the linear fit on logit scale.}
  \item{method}{either "gw.fit" or "logy.fit".  See the \code{method} argument.}
  \item{SS.null}{residual sum of squares from the null model.}
  \item{df.null}{residual degrees of freedom for the null model.}
  \item{SS.residual}{residual sum of squares from the fitted model.}
  
  In addition, non-null fits will have components \code{assign},
  \code{effects} and (unless not requested) \code{qr} relating to the
  linear fit, for use by extractor functions such as \code{summary}
  and \code{\link{effects}}.
}
\references{

Aitchison, J (1986).  \emph{The Statistical Analysis of Compositional Data}.
Chapman and Hall.  Republished and expanded as Aitchison, J. (2003).
\emph{The Statistical Analysis of Compositional Data}.  The Blackburn Press. 

Firth, D and Sammut, F (2023).  Analysis of composition on the original
scale of measurement.  \url{https://arxiv.org/abs/2312.10548}

Wedderburn, R. W. M. (1974). Quasi-likelihood functions, generalized
linear models, and the Gauss--Newton method.  \emph{Biometrika}
\bold{61}, 439--447.

}
\author{
David Firth (University of Warwick, UK) and Fiona Sammut (University of Malta)
}
\examples{
data(arctic_lake)
##
##  Fit the Aitchison-style multivariate linear model for 
##  dependence of logratios on log(depth):
sediments <- arctic_lake[, 1:3]
logdepth <- log(arctic_lake[, "depth"])
logratio_model <- colm(sediments ~ logdepth, ref = 3, method = "logratio_fit")
print(logratio_model)
## Those results replicate Aitchison (2003) [eqn (7.51) on page 165]. 
##
##  Fit the corresponding generalized Wedderburn logit 
##  model for dependence of composition on log(depth):
logit_model <- colm(sediments ~ logdepth, ref = 3)  ## method = "logit_fit"
print(logit_model)
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
