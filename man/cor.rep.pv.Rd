\name{cor.rep.pv}
\alias{cor.rep.pv}
\title{Multivariate correlation estimator from replicated gene microarray data}
\description{
This function computes p-values of correlation estimation from gene microarray data.
}
\usage{
cor.rep.pv(x, y=NULL, m, G)
}
\arguments{
  \item{x}{x is the gene microarray data in matrix OR data.frame format.}
  \item{y}{y optional. Used when x, y are vectors.}
  \item{m}{m is the number of replicates for each independent observations. Note: m*G == nrow(dat).}
  \item{G}{G is the number of the genes. Note: m*G == nrow(dat).}
}
\details{
This function computes p-values of the new multivariate correlation estimator from gene microarray data using permutation. It is used together with cor.rep(), and cor.rep.bootci() whenever necessary. 
}
\value{
This function returns a matrix of p-values corresponding to the matrix returned by cor.rep().
}
\references{
Zhu, D., Li, Y., and Hero, AO. Estimating gene expression correlation from replicated microarray data - A multivariate approach. {\it Submitted}
}
\author{Dongxiao Zhu (http://www-personal.umich.edu/~zhud)}

\seealso{\code{\link{cor.rep}}}
\examples{
# load GeneNT and GeneTS library
library(GeneNT)
d0 <- rnorm(100)
#sample size is set to 3 for quick demo purposes
for(l in 2:3)
d0 <- rbind(d0, rnorm(100))
d0<- t(d0)
M <- cor.rep(d0, m = 4, G= 25)
M.pv <- cor.rep.pv(d0, m = 4, G= 25)
}
\keyword{cluster}
\keyword{htest}
