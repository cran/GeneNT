cor.rep.bootci <- function(x, y=NULL, m, G, alpha)
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
Bootdist <- c(NA, 1000)
upperCI <- matrix(NA, D, D)
lowerCI <- matrix(NA, D, D)

for(i in 1:(D-1))
{
	for(j in (i+1):D)
	{
		subdat <- rbind(x[((i-1)*m + 1):((i-1)*m + m),], x[((j-1)*m + 1):((j-1)*m + m),])
		for(B in 1:1000)
		{
		subdatB <- subdat[, sample(seq(1:n), n, replace = TRUE)]
		muhatB <- c(rep(mean(subdatB[1:m,]), m), rep(mean(subdatB[(m+1):(2*m),]), m))
		ss0B <- (subdatB[,1] - muhatB)%o%(subdatB[,1] - muhatB)
		for (k in 2:n)
		ssB <- ss0B + (subdatB[,k] - muhatB)%o%(subdatB[,k] - muhatB)
		SigmahatB <- cov2cor(ssB/n)
		Bootdist[B] <- mean(SigmahatB[1:m,(m+1):(2*m)],SigmahatB[(m+1):(2*m),1:m])
		}
	upperCI[i,j] <- quantile(Bootdist, probs = c(alpha/2, 1-alpha/2))[1]
	lowerCI[i,j] <- quantile(Bootdist, probs = c(alpha/2, 1-alpha/2))[2]
	}

}
list(upperCI, lowerCI)
}
