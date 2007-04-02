\name{BEST.pearson}
\alias{BEST.pearson}
\title{BLAST type search for similar co-expressions with controlled FDR and MAS using Pearson correlation statistic}
\description{
This function implement the BLAST type two-stage screening procedure based on Pearson correlation coefficient. Specifying a pair of FDR and MAS criteria, and a seed gene, the algorithm screens a seeded gene cluster with controlled FDR and MAS. 
}
\usage{
BEST.pearson(gene.name, Q, cormin, method)
}
\arguments{
  \item{gene.name}{The "seed" gene name}
  \item{Q}{The significant level}
  \item{cormin}{The specified minimum acceptable strength of association measured by Pearson correlation statistic}
  \item{method}{The multiple comparison adjustment methods, eg, "BH" (for independent tests) and "BY" (for general dependent tests)}
}
\details{
The data matrix file must be in the right format. The first row (header) must be one shorter than the rest rows. The first column must be gene names.  
}
\value{
  The function returns a seeded gene cluster that satisfies the FDR and MAS criteria measured by Pearson correlation coefficient simultaneously. 
  \item{bpG1 }{The gene pairs that passes Stage I (FDR only) screening}
  \item{bpG2 }{The gene pairs that passed both Stage I (FDR) and II (MAS) screenings}
}
\references{
Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False 
Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{J Comput Biol}, \bold{12}, 1027-1043.

Fisher, R.A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32.
}
\author{Dongxiao Zhu (http://dongxiaozhu.googlepages.com/home)}

\seealso{\code{\link{BEST.kendall}}}
\examples{
## Load GeneNT library
library(GeneNT)

## load example dataset
data(dat) 
## Note, data matrix name has to be "dat"

## use (FDR, MAS) criteria (0.2, 0.5) and seed gene "GAL7" as example to screen gene pairs 
## use "BY" procedure to do multiple comparison adjustment
g4 <- BEST.pearson("GAL7", 0.2, 0.5, "BY")  
bpG1 <- g4$bpG1

## bpG2 contains gene pairs that passed two-stage screening
bpG2 <- g4$bpG2  
}
\keyword{cluster}
\keyword{htest}
