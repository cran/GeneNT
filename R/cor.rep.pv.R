cor.rep.pv <- function(x, y=NULL, m, G)
{

  if (is.data.frame(y)) 
        y <- as.matrix(y)
    else stopifnot(is.atomic(y))
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    else {
        stopifnot(is.atomic(x))
        if (!is.matrix(x)) {
            if (is.null(y)) 
                stop("supply both x and y or a matrix-like x")
            x <- as.vector(x)
        }
    }

p <- dim(x)[1]
n <- dim(x)[2]
D <- p/m
if(D != G)
stop("unqueal number of replicates or incorrect number of replicates or incorrect number of genes (random variables)")
PV <- matrix(NA, D, D)
nulldist <- rep(NA, factorial(n))
for(i in 1:(D-1))
{
	for(j in (i+1):D)
	{
	subdat <- rbind(x[((i-1)*m + 1):((i-1)*m + m),], x[((j-1)*m + 1):((j-1)*m + m),])
	muhat <- c(rep(mean(subdat[1:m,]), m), rep(mean(subdat[(m+1):(2*m),]), m))
	ss0 <- (subdat[,1] - muhat)%o%(subdat[,1] - muhat)
	for (k in 2:n)
	ss <- ss0 + (subdat[,k] - muhat)%o%(subdat[,k] - muhat)
	Sigmahat <- cov2cor(ss/n)
	cor0 <- mean(Sigmahat[1:m,(m+1):(2*m)],Sigmahat[(m+1):(2*m),1:m])		

		for(P in 1:(factorial(n)))
		{
		subdatP <- rbind(subdat[1:m, permutations(n)[P,]], subdat[(m+1):(2*m),])
		muhatP <- c(rep(mean(subdatP[1:m,]), m), rep(mean(subdatP[(m+1):(2*m),]), m))
		ss0P <- (subdatP[,1] - muhatP)%o%(subdatP[,1] - muhatP)
		for (k in 2:n)
		ssP <- ss0P + (subdatP[,k] - muhatP)%o%(subdatP[,k] - muhatP)
		SigmahatP <- cov2cor(ssP/n)
		nulldist[P] <- mean(SigmahatP[1:m,(m+1):(2*m)],SigmahatP[(m+1):(2*m),1:m])
		}
		PV[i,j] <- sum(abs(nulldist) > cor0)/factorial(n)
		PV[j,i] <- PV[i,j]
	}

}
diag(PV) <- 1
PV
}
