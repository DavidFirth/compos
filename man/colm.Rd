\name{cglm}
\alias{cglm}
\title{
Function to specify and fit multivariate regression models to compositional-response data
}
\description{
This function works in much the same way as \code{stats::glm}.  The models currently available are Aitchison-type
multivariate linear regression based on logratio transformation of the responses, and the 
Generalized Wedderburn logit model.
}
\usage{
cglm(formula, data, weights, subset, na.action = na.exclude, start = NULL, maxit = 100, 
   epsilon = 1e-08, trace = FALSE, method = "gw.fit", contrasts = NULL, ref = 1, 
   y.totals = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  a model formula.  See \code{\link{glm}} for more details.
}
  \item{data}{
  an optional data frame, list or environment as for \code{\link{glm}}.  If not found in \code{data}, 
  the variables for the model are taken from environment(formula), typically the environment from 
  which cglm is called.
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
  either \code{"gw.fit"} (for Generalized Wedderburn model) or \code{"logy.fit"} (for Aitchison-type multivariate
  linear model after logratio transformation).
  }
  \item{\dots}{
  arguments to pass to other functions.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A model object of class \Sexpr{"cglm"}, as a list with the following components:
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
  
  In addition, non-null fits will have components \code{assign}, \code{effects} and (unless not requested) \code{qr} relating to the linear fit, for use by extractor functions such as \code{summary} and 
  \code{\link{effects}}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
David Firth (University of Warwick, UK) and Fiona Sammut (University of Malta)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(ArcticLake)
##  Fit the Aitchison-style multivariate linear model for 
##  dependence of logratios on log(depth):
sediments <- ArcticLake[, 1:3]
logdepth <- log(ArcticLake[, "depth"])
cglm(sediments ~ logdepth, ref = 3, method = "logy.fit")
##  Fit the corresponding Generalized Wedderburn logit 
##  model for dependence of logratios on log(depth):
gw_model <- cglm(sediments ~ logdepth, ref = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
