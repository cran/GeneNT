\name{getBM}
\alias{getBM}
\title{Generate Pajek compatible matrix from screened gene pairs}
\description{
This function takes inputs of screened gene pairs based on BOTH Pearson correlation and Kendall correlation statistics, and generates Pajek compatible Boolean matrix that can be visulized using network visulization software such as Pajek. 
}
\usage{
getBM(pG2, kG2)
}
\arguments{
  \item{pG2}{Gene pairs that passed two-stage screening based on Pearson correlation coefficient.}
  \item{kG2}{Gene pairs that passed two-stage screening based on Kendall correlation coefficient.}
}
\details{
A Pajek compatible Boolean matrix "BMPajek.mat" and a R Boolean matrix "BM.tsv" will be exported to the working directory.
}
\value{
"BM.tsv" (R format Boolean matrix) and "BMPajek.mat" (Pajek format Boolean matrix) will be returned. 
}
\author{Dongxiao Zhu (\url{http://www-personal.umich.edu/~zhud})}

\seealso{\code{\link{row.names}}}
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
#g2 <- kendallfdrci(0.2, 0.5) #use (FDR, MAS) criteria (0.2, 0.5) as example to screen gene pairs 
#kG1 <- g2$kG1
#kG2 contains gene pairs that passed two-stage screening
#kG2 <- g2$kG2  
#generate Pajek compatible matrix to visualize network
#getBM(pG2, kG2) 
#clustering from network using network constraint clustering, for example. p=3.
#spclust(3, pG2, kG2) 
}
\keyword{cluster}
\keyword{htest}
