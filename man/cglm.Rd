\name{cglm}
\alias{cglm}
\title{
Function to specify and fit multivariate regression models to compositional-response data
}
\description{
This function works in much the same way as stats::glm.  The models currently available are Aitchison-type
multivariate linear regression based on logratio transformation of the responses, and the 
Generalized Wedderburn logit model.
}
\usage{
cglm(formula, data, weights, subset, na.action, start = NULL, maxit = 100, epsilon = 1e-08, trace = FALSE, model = TRUE, method = "gw.fit", y = TRUE, contrasts = NULL, ref = 1, y.totals = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  a model formula.  See stats::glm for more details.
}
  \item{data}{
  an optional data frame, list or environment as for stats::glm.  If not found in DATA, the variables for the
  model are taken from environmant(formula), typically the environment from which cglm is called.
}
  \item{weights}{
  an optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
}
  \item{subset}{
  an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{na.action}{
%%     ~~Describe \code{na.action} here~~
}
  \item{start}{
%%     ~~Describe \code{start} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{contrasts}{
%%     ~~Describe \code{contrasts} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
}
  \item{y.totals}{
%%     ~~Describe \code{y.totals} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(ArcticLake)
## Fit the Aitchison-style multivariate linear model for dependence of logratios on log(depth):
logdepth <- log(ArcticLake[, "depth"])
cglm(sediments ~ logdepth, ref = 3, method = "logy.fit")
## Fit the corresponding Generalized Wedderburn logit model for dependence of logratios on log(depth):
gw_model <- cglm(sediments ~ logdepth, ref = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
