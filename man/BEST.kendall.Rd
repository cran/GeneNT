\name{BEST.kendall}
\alias{BEST.kendall}
\title{BLAST type search for similar co-expressions with controlled FDR and MAS using Kendall correlation statistic}
\description{
This function implements the BLAST type two-stage screening procedure based on Kendall correlation statistic. Specifying a pair of FDR and MAS criteria, and a seed gene, the algorithm screens a seeded gene cluster with controlled FDR and MAS. 
}
\usage{
BEST.kendall(gene.name, Q, cormin)
}
\arguments{
  \item{gene.name}{The "seed" gene name}
  \item{Q}{The significant level}
  \item{cormin}{The specified minimum acceptable strength of association measured using Kendall correlation statistic}
}
\details{
The data matrix file must be in the right format. The first row (header) must be one shorter than the rest rows. The first column must be gene names.  
}
\value{
  The function returns a seeded gene cluster that satisfies the FDR and MAS constraints simultaneously using Kendall correlation statistic. 
  \item{bkG1}{The gene pairs that passes Stage I (FDR only) screening}
  \item{bkG2}{The gene pairs that passed both Stage I (FDR) and II (MAS) screenings}
}
\references{

Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{J of Comput Biol, in press}.

Hollander M. and Wolfe D.A. (1999). Nonparametric Statistical Methods, New York: Wiley.

}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{BEST.pearson}}}

\examples{
# load GeneNT and GeneTS library
library(GeneTS)
library(GeneNT)

#EITHER use the internal dataset
data(dat) 
#OR use the following if you want to import external data 
#dat <- read.table("gal.txt", h = T, row.names = 1) 
#Note, data matrix name has to be "dat"

#use (FDR, MAS) criteria (0.2, 0.5) and seed gene "GAL7" as example to screen gene pairs 
#g5 <- BEST.kendall("GAL7", 0.2, 0.5)  
#bkG1 <- g5$bkG1

#bkG2 contains gene pairs that passed the two-stage screening
#bkG2 <- g5$bkG2  
}
\keyword{cluster}
\keyword{htest}
