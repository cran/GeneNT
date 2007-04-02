\name{kendallfdrci}
\alias{kendallfdrci}
\title{Two-stage screening procedure based on Kendall correlation coefficient}
\description{
This function implement the two-stage screening procedure based on Kendall correlation coefficient. Specifying a pair of FDR and MAS criteria, the algorithm provides an initial co-expression discovery that controls only FDR, which is then followed by a second stage co-expression discovery which controls both FDR and MAS. 
}
\usage{
kendallfdrci(Q, cormin, method)
}
\arguments{
  \item{Q}{The significant level}
  \item{cormin}{The specified minimum acceptable strength of association measured using Kendall correlation coefficient}
  \item{method}{The multiple comparison adjustment methods, eg, "BH" (for independent tests) and "BY" (for general dependent tests)}
}
\details{
  The data matrix file must be in the right format. The first row must be one shorter 
  than the rest rows. The first column must be gene names.  
}
\value{
  The function returns a list of gene pairs that satisfies the FDR and MAS criteria simultaneously measured by Kendall correlation coefficient. 
  \item{kG1 }{The gene pairs that passes Stage I (FDR only) screening}
  \item{kG2 }{The gene pairs that passed both Stage I (FDR) and II (MAS) screenings}
}
\references{
Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False 
Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{J Comput Biol}, \bold{12}, 1027-1043.
Hollander M. and Wolfe D.A. (1999). Nonparametric Statistical Methods, New York: Wiley.
}
\author{Dongxiao Zhu (\url{http://dongxiaozhu.googlepages.com/home})}

\seealso{\code{\link{corfdrci}}}
\examples{
## load GeneNT library
library(GeneNT)

## load example dataset
data(dat) 
## Note, data matrix name has to be "dat"

g2 <- kendallfdrci(0.2, 0.5, "BY") 
kG1 <- g2$kG1
## kG2 contains gene pairs that passed two-stage screening
kG2 <- g2$kG2  
}
\keyword{cluster}
\keyword{htest}
