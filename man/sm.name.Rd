\name{sm.name}
\alias{sm.name}
\title{Extract gene names from data matrix}
\description{
This function extracts gene names from data matrix and use it for output.
}
\usage{
sm.name(m)
}
\arguments{
  \item{m}{M is a data matrix}
}
\value{
Returns a vector object of gene names.
}
\references{
Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{Submitted}.
}
\author{Dongxiao Zhu (http://www-personal.umich.edu/~zhud)}

\seealso{\code{\link{row.names}}}
\examples{
# load GeneNT and GeneTS library
library(GeneTS)
library(GeneNT)

#EITHER use the example dataset
data(dat) 
#OR use the following if you want to import external data 
#dat <- read.table("gal.txt", h = T, row.names = 1) 
#Note, data matrix name has to be "dat"
#calculate correlation matrix
M <- cor(t(dat)) 
#extract gene names (matrix object) from the correlation matrix
#NL <- sm.name(M) 
}
\keyword{cluster}
\keyword{htest}
