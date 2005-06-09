\name{BEST.pearson}
\alias{BEST.pearson}
\title{BLAST type search for similar co-expressions with controlled FDR and MAS using Pearson correlation statistic}
\description{
This function implement the BLAST type two-stage screening procedure based on Pearson correlation coefficient. Specifying a pair of FDR and MAS criteria, and a seed gene, the algorithm screens a seeded gene cluster with controlled FDR and MAS. 
}
\usage{
BEST.pearson(gene.name, Q, cormin, method = c("simple", "jackknife"))
}
\arguments{
  \item{gene.name}{The "seed" gene name}
  \item{Q}{The significant level}
  \item{cormin}{The specified minimum acceptable strength of association measured by Pearson correlation statistic}
  \item{method}{This argument can be set to either simple or jackknife if ONLY the correlation estimate is desired. It should be set to default (left blank) if p-values and CI's are desired.}
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
Zhu, D., Hero, A.O., Qin, Z.S. and Swaroop, A. High throughput screening of co-expressed gene pairs with controlled False Discovery Rate (FDR) and Minimum Acceptable Strength (MAS). \emph{J Comput Biol, in press}.

Fisher, R.A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small sample. \emph{Metron}, \bold{1}, 1--32.
}
\author{Dongxiao Zhu (http://www-personal.umich.edu/~zhud)}

\seealso{\code{\link{BEST.kendall}}}
\examples{
# load GeneNT and GeneTS library
library(GeneTS)
library(GeneNT)

#EITHER use the example dataset
data(dat) 
#OR use the following if you want to import external data 
#dat <- read.table("gal.txt", h = T, row.names = 1)  
#Note, data matrix name has to be "dat"

#use (FDR, MAS) criteria (0.2, 0.5) and seed gene "GAL7" as example to screen gene pairs 
#g4 <- BEST.pearson("GAL7", 0.2, 0.5)  
#bpG1 <- g4$bpG1

#bpG2 contains gene pairs that passed two-stage screening
#bpG2 <- g4$bpG2  
}
\keyword{cluster}
\keyword{htest}
