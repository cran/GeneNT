\name{priorclust}
\alias{priorclust}
\title{Regular "prior" clustering}
\description{
This function does regular clustering based on core R function hclust().
}
\usage{
priorclust(p)
}
\arguments{
  \item{p}{p is the scaling factor that need to be adjusted using biological prior knowledge.}
}
\details{
This function is written in comparison to the "posterior" clustering implemented as spclust(). 
}
\value{
This function returns a distance matrix that can be used by any distance based clustering software. 
}
\references{
Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{Submitted}.
}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{spclust}}}
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
#priorclust(3)
}
\keyword{cluster}
\keyword{htest}
