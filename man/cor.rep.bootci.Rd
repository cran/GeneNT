\name{cor.rep.bootci}
\alias{cor.rep.bootci}
\title{Bootstrap Confidence Interval for multivariate correlation estimator from replicated gene microarray data}
\description{
This function compute bootstrap confidence interval for correlation estimation from replicated gene microarray data.
}
\usage{
cor.rep.bootci(x, y=NULL, m, G, alpha)
}
\arguments{
  \item{x}{x is the gene microarray data in matrix OR data.frame format.}
  \item{y}{y (optional). Should be used when x, y are vectors.}
  \item{m}{m is the number of replicates for each independent observations. Note: m*G == nrow(dat).}
  \item{G}{G is the number of the genes. Note: m*G == nrow(dat).}
  \item{alpha}{alpha is the significance level.}
}
\details{
This function computes bootstrap confidence interval of the new multivariate correlation estimator at given significance level. It is used together with cor.rep(), and cor.rep.pv() whenever necessary. 
}
\value{
This function returns matrices of upper and lower bounds corresponding to the matrix returned by cor.rep().
}
\references{
Zhu, D., Li, Y., and Hero, AO. Estimating gene expression correlation from replicated microarray data - A multivariate approach. {\it Submitted}
}
\author{Dongxiao Zhu (http://www-personal.umich.edu/~zhud)}

\seealso{\code{\link{cor.rep}}}
\examples{
# load GeneNT and GeneTS library
library(GeneNT)
d0 <- rnorm(8)
#parameters are set to small values for quick demo purposes
for(l in 2:3)
d0 <- rbind(d0, rnorm(8))
d0<- t(d0)
M <- cor.rep(d0, m = 4, G= 2)
M.bootci <- cor.rep.bootci(d0, m = 4, G= 2, alpha = 0.05)
}
\keyword{cluster}
\keyword{htest}
