\name{kendall.confint}
\alias{kendall.confint}
\title{Asymptotic confidence intervals of Pearson correlation coefficient}
\description{
This function returns confidence intervals for a list of Kendall correlation coefficients given sample size N and significant level alpha.
}
\usage{
kendall.confint(x, y, alpha)
}
\arguments{
  \item{x}{i.i.d. observations of random variable X}
  \item{y}{i.i.d. observations of random variable Y}
  \item{alpha}{The significance level}
}
\details{
The function calculate approximate p-values and confidence intervals. See Hollander and Wolfe, 1999. Nonparametric Statistical Methods.
}
\value{
The function returns asymptotic lower and upper bound(s) of Kendall correlation coefficient(s). 
}
\references{
Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{Submitted}.

Hollander M. and Wolfe D.A. (1999). Nonparametric Statistical Methods, New York: Wiley.
}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{cor.confint}}}
\examples{
x <- runif(20, -1, 1)
y <- runif(20, -1, 1)
#kendall.confint(x,y,0.05)
}
\keyword{cluster}
\keyword{htest}
