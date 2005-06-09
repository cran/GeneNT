cor.rep <- function(x, y=NULL, m, G)
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
M <- matrix(NA, D, D)

for(i in 1:(D-1))
{
	for(j in (i+1):D)
	{
	subdat <- rbind(x[((i-1)*m + 1):((i-1)*m + m),], x[((j-1)*m + 1):((j-1)*m + m),])
	muhat <- c(rep(mean(subdat[1:m,]), m), rep(mean(subdat[(m+1):(2*m),]), m))
	
	#mudat <- c(rep(mean(dat[((i-1)*m + 1):((i-1)*m + m),]), m), rep(mean(dat[((j-1)*m + 1):((j-1)*m + m),]), m))           
	
	ss0 <- (subdat[,1] - muhat)%o%(subdat[,1] - muhat)
	for (k in 2:n)
	ss <- ss0 + (subdat[,k] - muhat)%o%(subdat[,k] - muhat)
	Sigmahat <- cov2cor(ss/n)
	M[i,j] <- mean(Sigmahat[1:m,(m+1):(2*m)],Sigmahat[(m+1):(2*m),1:m])
	M[j,i] <- M[i,j]
	}

}
diag(M) <- 1
M
}
