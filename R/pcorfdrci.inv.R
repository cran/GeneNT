pcorfdrci.inv <- function(pcormin)
{
	#dat <- read.table(file.name, h = T, row.names = 1)
        if (class(dat) != "matrix")
    	dat <- as.matrix(dat)		
	###step (1), get the p values
	gal.cor.m <- cor(t(dat))
	b.cor <- bagged.cor(gal.cor.m)
	b.pcor <- cor2pcor(b.cor)
	pcor.list <- sm2vec(b.pcor)
	kappa <- cor0.estimate.kappa(pcor.list) 
	pcov.p <- matrix(NA,nrow(dat), nrow(dat))

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
	p.index <- sm.indexes(pcov.p, diag = F)
	colnames(p.index) <- c("index1", "index2")
	p.list<- sm2vec(t(pcov.p), diag = F)
	G.results <- cbind(p.index, p.list, pcor.list)

	###step (2), sort the p values
	sort.index <- order(G.results[,3])
	G.results <- G.results[sort.index, ]
	
	###step (3),finding the min alpha  
	high <- array(NA, 99)
	low <- array(NA, 99)
	minalpha <- array(NA,length(pcor.list))
	fdrp <- array(NA,length(pcor.list))

	for (i in 1:length(pcor.list))
	{
	#within pcormins
	if(-pcormin <= G.results[i,4] && G.results[i,4] <= pcormin)
	minalpha[i] = 1

	#below neg pcormin
	if(G.results[i,4] < -pcormin)
	{
		for (al in 1:99)
		{
		g<- pcor.confint(G.results[i,4], kappa, al/100)
		high[al] <- g$conf.int2
		} 
		if ((max(high) >= -pcormin)&&(any(high < - pcormin)))
		minalpha[i] <- min(which(low > pcormin))/100

		if(max(high) < -pcormin)
		minalpha[i] <- 0

		if(!any(high < - pcormin))
		minalpha[i] <- 1
	}

	#above pos pcormin
	if(G.results[i,4] > pcormin)
	{
		for (al in 1:99)
		{
		g<- pcor.confint(G.results[i,4], kappa, al/100)
		low[al] <- g$conf.int1
		} 
		if((min(low) <= pcormin)&&(any(low > pcormin)))
		minalpha[i] <- min(which(low > pcormin))/100

		if((min(low) > pcormin))
		minalpha[i] <- 0

		if(!any(low > pcormin))
		minalpha[i] <- 1
	}
	}  
	###step (4), compute the index
	boo <- array(NA, length(pcor.list))
	fdrp <- array(NA, length(pcor.list)) 
	for (j in 1:length(pcor.list))
	{
		for (k in 1:length(pcor.list))
		{
		boo[k] <- (G.results[k,3]*k/length(pcor.list) <= minalpha[j])
		}
		fdrp[j] <- ifelse(sum(boo) == 0, G.results[1,3], G.results[sum(boo),3])
	}
        list(fdrp = fdrp)
}

