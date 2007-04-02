kendallfdrci <- function(Q, cormin, method)
{
    if (class(dat) != "matrix")
    dat <- as.matrix(dat)		
	  gal.cor.m <- cor(t(dat), method = "kendall")
	  cov.p <- matrix(NA, nrow(dat), nrow(dat))
	  for (i in 1:(nrow(dat)-1))
	  {
		    for (j in (i+1):nrow(dat))
		    {
		    g<- cor.test(dat[i,], dat[j,], alternative = "two.sided", method = "kendall", exact = FALSE)
		    cov.p[i,j] <- g$p.value
		    cov.p[j,i] <- cov.p[i,j]
		    }
	  }
	  ###Calculate p-values for all pairs
	  p.name <- sm.name(gal.cor.m)
	  colnames(p.name) <- c("gene1", "gene2")
	  p.name <- data.frame(p.name)
    p.list <- as.vector(cov.p[lower.tri(cov.p, diag = FALSE)])
    cor.list <- as.vector(gal.cor.m[lower.tri(gal.cor.m, diag = FALSE)]) 
	  p.results <- cbind(p.name, cor.list, p.list)

	  ###estimate fdr using either BH or BY procedures
    adjp.list <- p.adjust(p.list, method=method)
    adjp.results <- cbind(p.results, adjp.list)

	  ###use stepdown procedure to generate G1 set at certain FDR
	  G1.results <- adjp.results[(adjp.results[,5] < Q),]
	  sort.idx <- order(-abs(G1.results[, 3]))
	  G1.results <- G1.results[sort.idx, ]
	  G1 <- nrow(G1.results)
	  pq <- G1/length(p.list)

	  ###calculate confidence intervals for G1
	  lower <- rep(NA, G1)
	  higher <- rep(NA, G1)
	  alpha = Q*pq
	
	  if(G1 == 0)
    warning("Stage I screening returns no results. The Q value may be set to conservative (low)!")		
	  for(i in 1:G1)
	  {
        x <- which(row.names(dat) == G1.results[i,1])
	      y <- which(row.names(dat) == G1.results[i,2])
	      g <- kendall.confint(dat[x,], dat[y,], alpha)
	      lower[i] <- g[1]
	      ifelse(g[2] >1, higher[i] <- 1, higher[i] <- g[2])
	  }
  	kG1 <- cbind(G1.results, lower, higher) 
    indxg2 <- apply(kG1, 1, function(x) ifelse((x[6] > cormin) || (x[7] < -cormin), indxg2<- TRUE, indxg2 <- FALSE))
    kG2 <- kG1[indxg2,]
  	write.table(kG2, sep = "\t", file = "kG2.txt")
  	write.table(kG1, sep = "\t", file = "kG1.txt")
  	cat("The screened pairs are now in your working directory. ")
  	list(kG1=kG1, kG2=kG2) 
}

