\name{pcorfdrci.inv}
\alias{pcorfdrci.inv}
\title{Inverse screening procedure based on Kendall correlation coefficient}
\description{
The inverse screening procedure is presented in our the paper below. 
It allows the FDR p-value of a gene pair's strength of association measured 
using partial correlation coefficient to be computed
}
\usage{
pcorfdrci.inv(pcormin)
}
\arguments{
  \item{pcormin}{partial correlation threshold}
}
\details{
Please refer to Zhu et al. below for detail. 
}
\value{
This function returns a list of fdr adjusted p-values.
}
\references{
Fisher, R.A. (1921). On the 'probable error' of a coefficient of correlation deduced from
  a small sample. \emph{Metron}, \bold{1}, 1--32.

Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed 
gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{Submitted}.

Schäfer, J., and K. Strimmer. (2004) An empirical Bayes approach to inferring 
large-scale gene association networks. \emph{Bioinformatics}, \bold{1}, 1-13.
}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{corfdrci.inv}}}
\examples{
# load GeneNT and GeneTS library
library(GeneTS)
library(GeneNT)

#EITHER use the example dataset
data(dat) 
#OR use the following if you want to import external data 
#dat <- read.table("gal.txt", h = T, row.names = 1) 
#Note, data matrix name has to be "dat"
#Subset data to run faster
dat <- dat[1:10,]
#Set MAS level pcormin = 0.8, calculate fdr adjusted p-values. 
#fdrp8 <- pcorfdrci.inv(0.8)                                   
}
\keyword{cluster}
\keyword{htest}
