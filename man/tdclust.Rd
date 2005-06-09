\name{tdclust}
\alias{tdclust}
\title{Traditional clustering}
\description{
This function does traditional clustering based on core R function hclust().
}
\usage{
tdclust(p)
}
\arguments{
  \item{p}{p is the exponential tuning factor. The default value is 1, and it can be set to other positive integers when necessary.}
}
\details{
This function is written in comparison to the network constrained clustering implemented as ncclust(). 
}
\value{
This function returns a distance matrix that can be used by any distance based clustering software. 
}
\references{
Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{J. Comput. Biol, in press}.
}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{ncclust}}}
\examples{
# load GeneNT and GeneTS library
library(GeneTS)
library(GeneNT)

#EITHER use the example dataset
data(dat) 
#OR use the following if you want to import external data 
#dat <- read.table("gal.txt", h = T, row.names = 1) 
#Note, data matrix name has to be "dat"
#Regular hierarchical clustering, for example, p = 3.
#tdclust(3)
}
\keyword{cluster}
\keyword{htest}
