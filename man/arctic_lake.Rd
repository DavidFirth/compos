\name{arctic_lake}
\alias{arctic_lake}
\docType{data}
\title{
  Stanwell-Fletcher Lake Sediment Data
}
\description{
  The depth and composition --- into (sand, silt, clay) --- of 39 sediment samples are
  recorded in a 39 by 4 numeric matrix.  The compositions all total 100
  to within rounding error.  Rows of the matrix are ordered and named as
  in Aitchison (1986, 2003); they are in order of depth apart from row 18.  Two
  data-transcription errors that appeared in Aitchison (1986, 2003) have
  been corrected here. 
}
\usage{data("arctic_lake")}
\format{
  A 39 by 4 numeric matrix.
}
\source{
  Coakley (1966, Table 5).

  The data were later reproduced in Coakley and Rust
(1968, Table 1) but with one transcription error: the "sand" part of
sediment "S4" was erroneously recorded there as 52.2 rather than 52.5.
The data were also reproduced in Aitchison (1986, 2003) as "Data 5", where the
error from Coakley and Rust (1968) was repeated.  A further
transcription error was introduced in Aitchison (1986, 2003): the
"silt" part of sediment "S24" was erroneously recorded there as 54.7 rather than 54.2.
Both of these errors have been corrected here.
}
\references{

Aitchison, J (1986).  \emph{The Statistical Analysis of Compositional Data}.
Chapman and Hall.

Aitchison, J. (2003). \emph{The Statistical Analysis of Compositional Data}.  The Blackburn Press. 

Coakley, J P (1966).  History and bottom sediments of Stanwell-Fletcher
Lake, Somerset Island, N.W.T.  MSc Thesis, University of Ottawa. \url{https://ruor.uottawa.ca/bitstream/10393/6966/1/MK01764.PDF}

Coakley, J P and Rust, B R (1968). Sedimentation in an Arctic lake.
\emph{Journal of Sedimentary Petrology} \bold{38}(4), 1290--1300. \doi{10.1306/74D71B59-2B21-11D7-8648000102C1865D}

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
##  Those results essentially replicate Aitchison (2003) [eqn (7.51) on page 165].
##
##  Fit the corresponding 'generalized Wedderburn' logit model for 
##  dependence of composition on log(depth):
logit_model <- colm(sediments ~ logdepth, ref = 3, method = "logit_fit")
print(logit_model)
}
\keyword{datasets}
