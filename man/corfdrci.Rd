\name{corfdrci}
\alias{corfdrci}
\title{Two-stage screening procedure for screening gene pairs based on Pearson correlation statistic}
\description{
This function implements the two-stage screening procedure based on Pearson correlation coefficient to screen similarly co-expressed gene pairs. Specifying a pair of FDR and MAS criteria, the algorithm provides an initial co-expression discovery that controls only FDR, which is then followed by a second stage co-expression discovery which controls both FDR and MAS.}
\usage{
corfdrci(Q, cormin)
}
\arguments{
  \item{Q}{The significant level}
  \item{cormin}{The specified minimum acceptable strength of association measured using Pearson correlation coefficient}
}
\details{
  The data matrix file must be in the right format. The first row (header) must be one shorter than the rest rows. The first column must be gene names.  
}
\value{
  The function returns a list of gene pairs that satisfies the FDR and MAS criteria simultaneously measured by Pearson correlation statistic. 
  \item{pG1 }{The gene pairs that passes Stage I (FDR only) screening}
  \item{pG2 }{The gene pairs that passed both Stage I (FDR) and II (MAS) screenings}
}
\references{

Fisher, R.A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32.

Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{Submitted}.

Schäfer, J., and K. Strimmer. (2004) An empirical Bayes approach to inferring large-scale gene association networks. \emph{Bioinformatics}, \bold{1}, 1-13.
}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{kendallfdrci}}}
\examples{

# load GeneNT and GeneTS library
library(GeneTS)
library(GeneNT)
library(e1071)

#EITHER use the example dataset
data(dat) 
#OR use the following if you want to import external data 
#dat <- read.table("gal.txt", h = T, row.names = 1) 
#Note, data matrix name has to be "dat"

#use (FDR, MAS) criteria (0.2, 0.5) as example to screen gene pairs 
#g1 <- corfdrci(0.2, 0.5) 
#pG1 <- g1$pG1
#pG2 contains gene pairs that passed two-stage screening
#pG2 <- g1$pG2  
#use (FDR, MAS) criteria (0.2, 0.5) as example to screen gene pairs 
#g2 <- kendallfdrci(0.2, 0.5) 
#kG1 <- g2$kG1
#kG2 contains gene pairs that passed two-stage screening
#kG2 <- g2$kG2  
#generate Pajek compatible matrix to visualize network
#getBM(pG2, kG2) 
#clustering from network using network constraint clustering, for example, p = 3.
#spclust(3, pG2, kG2) 
}
\keyword{cluster}
\keyword{htest}
