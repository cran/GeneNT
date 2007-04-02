\name{cor.confint}
\alias{cor.confint}
\title{Asymptotic confidence intervals of Pearson correlation coefficient}
\description{
This function returns confidence intervals for a list of Pearson correlation estimates given a sample size N and significant level alpha.}
\usage{
cor.confint(cor, N, alpha)
}
\arguments{
  \item{cor}{Pearson correlation coefficient estimate, can be single value or vector.}
  \item{N}{Sample size.}
  \item{alpha}{The significant level}
}
\details{
The function use Fisher transformation to approximate p-values and confidence intervals.
}
\value{
  The function returns asymptotic lower and upper bound(s) of Pearson correlation coefficient(s). 
  \item{conf.int1 }{Lower bound(s) of Pearson correlation coefficient(s)}
  \item{conf.int2 }{Upper bound(s) of Pearson correlation coefficient(s)}
}
\references{

Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{J Comput Biol}.

Schäfer, J., and K. Strimmer. (2004) An empirical Bayes approach to inferring large-scale gene association networks. \emph{Bioinformatics}, \bold{1}, 1-13.

}
\author{Dongxiao Zhu (\url{http://dongxiaozhu.googlepages.com/home})}

\seealso{\code{\link{kendall.confint}}}
\examples{
## simulate a vector of correlation coefficients 
cor <- runif(20, -1, 1) 
g6 <- cor.confint(cor, 20, 0.05)
g6$conf.int1
g6$conf.int2
}
\keyword{cluster}
\keyword{htest}
