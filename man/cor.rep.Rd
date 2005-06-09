\name{cor.rep}
\alias{cor.rep}
\title{Estimating gene expression correlation from replicated microarray data}
\description{
This function estimates the gene expression correlation from replicated microarray data. See exmaple for data format. 
}
\usage{
cor.rep(x, y=NULL, m, G)
}
\arguments{
  \item{x}{x is the gene microarray data in matrix OR data.frame format.}
  \item{y}{y Required only both x,y are vectors.}
  \item{m}{m is the number of replicates for each independent observations. Note: m*G == nrow(dat).}
  \item{G}{G is the number of the genes. Note: m*G == nrow(dat).}
}
\details{
This function implements the new multivariate correlation estimator from gene microarray data. The estimator takes
all observations rather than averaging to gain better accuracy especially for noisy microarray data. See the manuscript
for detail. 
}
\value{
This function returns a correlation matrix that corresponds to R function cor().
}
\references{
Zhu, D., Li, Y., and Hero, AO. Estimating gene expression correlation from replicated microarray data - A multivariate approach. {\it Submitted}
}
\author{Dongxiao Zhu (http://www-personal.umich.edu/~zhud)}

\seealso{\code{\link{cor}}}
\examples{
# load GeneNT and GeneTS library
library(GeneNT)
d0 <- rnorm(100)
for(l in 2:10)
d0 <- rbind(d0, rnorm(100))
d0<- t(d0)
M <- cor.rep(d0, m = 4, G= 25)
}
\keyword{cluster}
\keyword{htest}
