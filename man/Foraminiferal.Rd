\name{Foraminiferal}
\alias{Foraminiferal}
\docType{data}
\title{
Foraminiferal compositional data from Aitchison (1986)
}
\description{
The Foraminiferal data is made up of 30 compositions of foraminifer, a single-cell marine micro-organism, with the compositions being recorded at different water depths (in metres) with depth varying from 1m to 30m.  The compositional variables considered in this dataset are Neogloboquadrina atlantica (Na), Neogloboquadrina pachyderma (Np), Globorotalia obesa (Go) and Globigerinoides triloba (Gt).  
}
\usage{data("Foraminiferal")}
\format{
  A 30 by 5 numeric matrix.
}
\source{
Aitchison (1986), dataset 34. 
}
\references{
Aitchison, J (1986).  \emph{The Statistical Analysis of Compositional Data}.  Chapman and Hall.
}
\examples{
data(Foraminiferal)
head(Foraminiferal)
foraminifer <- Foraminiferal[, c("Na", "Np", "Go", "Gt")]
rowSums(foraminifer)
depth  <- Foraminiferal[, "depth"]
##  Fit the Generalized Wedderburn logit model for dependence of composition on depth:
logit_model <- colm(foraminifer ~ depth, ref = 4)
}
\keyword{datasets}
