corfdrci <- function(Q, cormin)
{
	if (class(dat) != "matrix")
    	dat <- as.matrix(dat)		
	gal.cor.m <- cor(t(dat))
	cov.p <- matrix(NA,nrow(dat), nrow(dat))
	for (i in 1:(nrow(dat)-1))
	{
		for (j in (i+1):nrow(dat))
		{
		g<- cor.test(dat[i,], dat[j,], alternative = "two.sided", method = "pearson")
		cov.p[i,j] <- g$p.value
		}
	}
	###Calculate p-values for all pairs
	p.name <- sm.name(gal.cor.m)
	colnames(p.name) <- c("gene1", "gene2")
	p.name <- data.frame(p.name)
	p.index <- sm.indexes(cov.p, diag = F)
	colnames(p.index) <- c("index1", "index2")
	p.list<- sm2vec(t(cov.p), diag = F)
	cor.list <- sm2vec(t(gal.cor.m), diag = F) 
	p.results <- cbind(p.index, p.name, cor.list, p.list)
	
	###calculate q values
	fdr.out <- fdr.control(p.list, Q)
	q.list <- fdr.out$qvalues
	q.results <- cbind(p.results, q.list)

	###use stepdown procedure to generate G1 set at certain FDR
	G1.results <- q.results[(q.results[,7] < Q),]
	sort.idx <- order(-abs(G1.results[, 5]))
	G1.results <- G1.results[sort.idx, ]
	G1 <- nrow(G1.results)
	pq <- G1/length(p.list)

	###calculate confidence intervals for G1
	lower <- array(NA, G1)
	higher <- array(NA, G1)
	alpha = Q*pq

	if(G1 == 0)
	warning("Stage I screening returns no results. The Q value may be set to conservative (low)!")		
	for(i in 1:G1)
	{
	x <- G1.results[i,1]
	y <- G1.results[i,2]
	r<- cor(dat[x,], dat[y,])
	g <- cor.confint(r, ncol(dat), alpha)
	lower[i] <- g$conf.int1
	higher[i] <- g$conf.int2
	}
	lower <- as.matrix(lower)
	lower <- data.frame(lower)
	higher <- as.matrix(higher)
	higher <- data.frame(higher)
	pG1 <- cbind(G1.results, lower, higher) 
	indxg2 <- apply(pG1,1,function(x) ifelse(((as.numeric(x[8]) > cormin) | (as.numeric(x[9]) < -cormin)), indxg2<- T, indxg2 <- F))
	pG2 <- pG1[indxg2,]
	write.table(pG2, sep = "\t", file = "pG2.txt")
	write.table(pG1, sep = "\t", file = "pG1.txt")
        cat("The screened pairs are now in your working directory. ")
	list(pG1=pG1, pG2=pG2)
}

