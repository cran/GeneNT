\name{pcorfdrci}
\alias{pcorfdrci}
\title{Two-stage screening procedure based on partial correlation coefficient in Graphic Gaussian Model framework}
\description{
This function implement the two-stage screening procedure based on partial correlation coefficient in the framework of Gaussian Graphic Model. Specifying a pair of FDR and MAS criteria, the algorithm provides an initial co-expression discovery that controls only FDR, which is then followed by a second stage co-expression discovery which controls both FDR and MAS. 
}
\usage{
pcorfdrci(Q, pcormin)
}
\arguments{
  \item{Q}{The significant level}
  \item{pcormin}{The specified minimum acceptable strength of association measured using partial correlation coefficient}
}
\details{
  The data matrix file must be in the right format. The first row (header) must be one shorter than the rest rows. The first column must be gene names.  
}
\value{
  The function returns a list of gene pairs that satisfies the FDR and MAS criteria simultaneously measured by partial correlation coefficient. 
  \item{pG1 }{The gene pairs that passes Stage I (FDR only) screening}
  \item{pG2 }{The gene pairs that passed both Stage I (FDR) and II (MAS) screenings}
}
\references{
Fisher, R.A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32.

Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{Submitted}.

Schäfer, J., and K. Strimmer. (2004) An empirical Bayes approach to inferring large-scale gene association networks. \emph{Bioinformatics}, \bold{1}, 1-13.
}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{corfdrci}}}
\examples{
# load GeneNT and GeneTS library
library(GeneTS)
library(GeneNT)

#EITHER use the example dataset
data(dat) 
#OR use the following if you want to import external data 
#dat <- read.table("gal.txt", h = T, row.names = 1) 
#Note, data matrix name has to be "dat"

#use (FDR, MAS) criteria (0.2, 0.2) as example to screen gene pairs 
#g3 <- pcorfdrci(0.2, 0.2) 
#G1 <- g3$G1.all
#G2 is the dataset containing gene pairs that passed two-stage screening
#G2 <- g3$G2  
}
\keyword{cluster}
\keyword{htest}
