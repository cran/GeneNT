pcorfdrci <- function(Q, pcormin)
{
	#dat <- read.table(file.name, h = T, row.names = 1)
        if (class(dat) != "matrix")
    	dat <- as.matrix(dat)	
	gal.cor.m <- cor(t(dat))
	b.pcor <- ggm.estimate.pcor(t(dat), method = "bagged.pcor")
	pcor.list <- sm2vec(b.pcor)
	
	pcov.p <- matrix(NA,nrow(dat), nrow(dat))
	kappa <- cor.fit.mixture(pcor.list)$kappa 
	for (i in 1:(nrow(dat)-1))
	{
		for (j in (i+1):nrow(dat))
		{
		r<- b.pcor[i,j]
		pval<- cor0.test(r, kappa, method = "ztransform")
		pcov.p[i,j] <- pval
		}
	}	
	###Calculate p values for all pairs
	p.name <- sm.name(gal.cor.m)
	colnames(p.name) <- c("gene1", "gene2")
	p.name <- data.frame(p.name)
	p.index <- sm.indexes(pcov.p, diag = F)
	colnames(p.index) <- c("index1", "index2")
	p.list<- sm2vec(t(pcov.p), diag = F)
	p.results <- cbind(p.index,p.name, pcor.list, p.list)

	###calculate q values
	fdr.out <- fdr.control(p.list, Q)
	q.list <- fdr.out$qvalues
	q.results <- cbind(p.results, q.list)

	###stepdown procedure to generate G1 subset controling FDR
	G1.results <- q.results[q.results[,7] < Q,]
	sort.idx <- order(-abs(G1.results[, 5]))
	G1.results <- G1.results[sort.idx, ]
	G1 <- nrow(G1.results)
	pq <- G1/length(p.list)

	###calculate confidence intervals for G1
	lower <- array(NA, G1)
	higher <- array(NA, G1)
	alpha = Q*pq

	for(i in 1:G1)
	{
		x <- G1.results[i,1]
		y <- G1.results[i,2]
		pcor <- b.pcor[x,y]
		g<- pcor.confint(pcor, kappa, alpha)

		lower[i] <- g$conf.int1
		higher[i] <- g$conf.int2
	}
	
	lower <- as.matrix(lower)
	lower <- data.frame(lower)
	higher <- as.matrix(higher)
	higher <- data.frame(higher)

	G1.all <- cbind(G1.results, lower, higher) 
	indxg2 <- apply(G1.all,1,function(x) ifelse(((as.numeric(x[8]) > pcormin) | (as.numeric(x[9]) < - pcormin)), indxg2<- T, indxg2 <- F))
	G2 <- G1.all[indxg2,]
	write.table(G2, sep = "\t", file = "G2.txt")
	write.table(G1.all, sep = "\t", file = "G1.txt")
	cat("The screened pairs are now in your working directory. ")
}

