\name{pcor.confint}
\alias{pcor.confint}
\title{Asymptotic confidence intervals of partial correlation coefficient}
\description{
This function returns confidence intervals for a list of partial
correlation coefficients given empirical sample size kappa and significant level alpha.}
}
\usage{
pcor.confint(pcor, kappa, alpha)
}
\arguments{
  \item{pcor}{Partial correlation coefficient(s)}
  \item{kappa}{Empirical sample size}
  \item{alpha}{The significant level}
}
\details{
For details about the empirical sample size, refer to Schäfer, J., and K. Strimmer, 2004, Bioinformatics, 1:1-13.
Briefly, it refers sample size that are estimated from data using a two-component 
mixture model rather than from sample size in the situation of "small N, large p" in 
the Gaussian Graphic Model framework. 
}
\value{
  The function returns asymptotic lower and upper bound(s) of Pearson correlation coefficient(s). 
  \item{conf.int1 }{Lower bound(s) of Pearson correlation coefficient(s)}
  \item{conf.int2 }{Upper bound(s) of Pearson correlation coefficient(s)}
}
\references{

Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of 
co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum 
Acceptable Strength (MAS). \emph{Submitted}.

Schäfer, J., and K. Strimmer. (2004) An empirical Bayes approach to inferring 
large-scale gene association networks. \emph{Bioinformatics}, \bold{1}, 1-13.

}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{cor.confint}}}
\examples{
pcor <- runif(20, -1, 1) #simulate a vector of correlation coefficients 
#g7 <- pcor.confint(pcor, 20, 0.05)
#g7$conf.int1
#g7$conf.int2
}
\keyword{cluster}
\keyword{htest}
