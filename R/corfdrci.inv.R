corfdrci.inv <- function(cormin)
{
	#dat <- read.table(file.name, h = T, row.names = 1)
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
	###Calculate p values for all pairs
	p.index <- sm.indexes(cov.p, diag = F)
	colnames(p.index) <- c("gene1", "gene2")
	cor.list <- sm2vec(t(gal.cor.m), diag = F) 
	p.list<- sm2vec(t(cov.p), diag = F)
	G.results <- cbind(p.index, p.list, cor.list)

	###step (2), sort the p values
	sort.index <- order(G.results[,3])
	G.results <- G.results[sort.index, ]
	###step (3),finding the min alpha  
	high <- array(NA, 99)
	low <- array(NA, 99)
	minalpha <- array(NA,length(cor.list))
	fdrp <- array(NA,length(cor.list))

	for (i in 1:length(cor.list))
	{
		x <- G.results[i,1]
		y <- G.results[i,2]

		#within cormin
		if(-cormin <= G.results[i,4] && G.results[i,4] <= cormin)
		minalpha[i] = 1
		
		#below neg cormin
		if(G.results[i,4] < -cormin)
		{
		for (al in 1:99)
		high[al] <- (cor.test(dat[x,], dat[y,], alternative = "two.sided", method = "pearson", conf.level = al/100))$conf.int[2] 
		
		if ((max(high) >= -cormin)&&(any(high < - cormin)))
		minalpha[i] <- min(which(high < -cormin))/100
		
		if(max(high) < -cormin)
		minalpha[i] <- 0
		
		if(!any(high < -cormin))
		minalpha[i] <- 1
		}

		#above pos cormin
		if(G.results[i,4] > cormin)
		{
		for (al in 1:99)
		low[al] <- (cor.test(dat[x,], dat[y,], alternative = "two.sided", method = "pearson", conf.level = al/100))$conf.int[1] 

		if((min(low) <= cormin)&&(any(low > cormin)))
		minalpha[i] <- (100- max(which(low > cormin)))/100
		
		if((min(low) > cormin))
		minalpha[i] <- 0

		if(!any(low > cormin))
		minalpha[i] <- 1
		}
	}  
	###step (4), compute the index
	boo <- array(NA, length(cor.list))
	fdrp <- array(NA, length(cor.list)) 
	for (j in 1:length(cor.list))
	{
		for (k in 1:length(cor.list))
		{
		boo[k] <- (G.results[k,3]*k/length(cor.list) <= minalpha[j])
		}
	fdrp[j] <- ifelse(sum(boo) == 0, G.results[1,3], G.results[sum(boo),3])
	}
	list(fdrp = fdrp)
}

