\name{dat}
\alias{dat}
\docType{data}
\title{Part of galactose metabolism data. Ideker et al. Science, 292:929-934}
\description{
The dataset contains top 100 most differentially expressed genes. 
}
\usage{data(dat)}
\details{
The whole dataset can be downloaded from supplemental website of the paper.
}
\source{
Science, 292:929-934, 2001.
}
\references{
Ideker T, Thorsson V, Ranish JA, Christmas R, Buhler J, Eng JK, Bumgarner R, Goodlett DR, Aebersold R, Hood L.(2001) Integrated genomic and proteomic analyses of a systematically perturbed metabolic network. \emph{Science}. \bold{292}, 929-934.
}
\examples{
data(dat)
boxplot(dat, xlab = "different experimental conditions", ylab = "log2 ratio of cy3/cy5 intensities") 
}
\keyword{datasets}
