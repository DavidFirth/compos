\name{ArcticLake}
\alias{ArcticLake}
\docType{data}
\title{
  Arctic lake sediment data from Aitchison (1986)
}
\description{
  The depth and composition --- into (sand, silt, clay) --- of 39 sediment samples are
  recorded in a 39 by 4 numeric matrix.  The compositions all total approximately 100
  percent.
}
\usage{data("ArcticLake")}
\format{
  A 39 by 4 numeric matrix.
}
\source{
Aitchison (1986), dataset 5.  Originally from a paper of Coakley and Rust, Journal of Sedimentary Petrology, 1968.
}
\references{
Aitchison, J (1986).  The Statistical Analysis of Compositional data.  Chapman and Hall.
}
\examples{
data(ArcticLake)
head(ArcticLake)
sediments <- ArcticLake[, c("sand", "silt", "clay")]
rowSums(sediments)
##  Fit the Aitchison-style multivariate linear model for dependence 
##  of logratios on log(depth):
logdepth <- log(ArcticLake[, "depth"])
cglm(sediments ~ logdepth, ref = 3, method = "logy.fit")
##  Fit the corresponding Generalized Wedderburn logit model for 
##  dependence of logratios on log(depth):
gw_model <- cglm(sediments ~ logdepth, ref = 3)
}
\keyword{datasets}
