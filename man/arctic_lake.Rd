\name{arctic_lake}
\alias{arctic_lake}
\docType{data}
\title{
  Arctic lake sediment data from Aitchison (1986)
}
\description{
  The depth and composition --- into (sand, silt, clay) --- of 39 sediment samples are
  recorded in a 39 by 4 numeric matrix.  The compositions all total approximately 100
  percent.
}
\usage{data("arctic_lake")}
\format{
  A 39 by 4 numeric matrix.
}
\source{
Aitchison (1986), dataset 5.  Originally from a paper of Coakley and Rust, Journal of Sedimentary Petrology, 1968.
}
\references{

Aitchison, J (1986).  \emph{The Statistical Analysis of Compositional Data}.
Chapman and Hall.  Republished and expanded as Aitchison, J. (2003).
\emph{The Statistical Analysis of Compositional Data}.  The Blackburn Press. 

Coakley, J P and Rust, B R (1968). Sedimentation in an Arctic lake.
\emph{Journal of Sedimentary Petrology} \bold{38}(4), 1290--1300.

}
\examples{
data(arctic_lake)
head(arctic_lake)
sediments <- arctic_lake[, c("sand", "silt", "clay")]
rowSums(sediments)
##
##  Fit the Aitchison-style multivariate linear model for dependence 
##  of logratios on log(depth):
logdepth <- log(arctic_lake[, "depth"])
logratio_model <- colm(sediments ~ logdepth, ref = 3, method = "logratio_fit")
print(logratio_model)
## Those results replicate Aitchison (2003) [eqn (7.51) on page 165]. 
##
##  Fit the corresponding Generalized Wedderburn logit model for 
##  dependence of composition on log(depth):
logit_model <- colm(sediments ~ logdepth, ref = 3)
print(logit_model)
}
\keyword{datasets}
