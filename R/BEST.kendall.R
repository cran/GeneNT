BEST.kendall <- function(gene.name, Q, cormin, method)
{
    if (class(dat) != "matrix")
    dat <- as.matrix(dat)    
    qy <- which(row.names(dat)==gene.name)
    p.list <- array(NA, nrow(dat))
    cor.list <- array(NA, nrow(dat))
    for (j in 1:nrow(dat))
    {
        g <- cor.test(dat[qy,], dat[j,], alternative = "two.sided", method = "kendall", exact = FALSE)
        cor.list[j] <- cor(dat[qy,], dat[j,], method = "kendall")
        p.list[j] <- g$p.value
    }
    ###Calculate p-values for all pairs
    p.name1 <- rep(row.names(dat)[qy], length(cor.list))
    p.name2 <- row.names(dat)
    p.name <- cbind(p.name1, p.name2)
    colnames(p.name) <- c("gene1", "gene2")
    p.name <- data.frame(p.name)
    p.results <- cbind(p.name, cor.list, p.list)
    ###calculate adjusted pvalues using BH or BY procedures
    adjp.list <- p.adjust(p.list, method=method)
    adjp.results <- cbind(p.results, adjp.list)
    ###use stepdown procedure to generate G1 set controlling FDR
    G1.results <- adjp.results[(adjp.results[,5] < Q),]
    sort.idx <- order(-abs(G1.results[, 3]))
    G1.results <- G1.results[sort.idx, ]
    G1 <- nrow(G1.results)
    pq <- G1/length(p.list)
    alpha = Q*pq
    ###calculate confidence intervals for G1
    lower <- array(NA, G1)
    higher <- array(NA, G1)
    
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
    bkG1 <- cbind(G1.results, lower, higher) 
    indxg2 <- apply(bkG1,1,function(x) ifelse((x[6] > cormin) || (x[7] < -cormin), indxg2<- TRUE, indxg2 <- FALSE))
    bkG2 <- bkG1[indxg2,]
    write.table(bkG2, sep = "\t", file = "bkG2.txt")
    write.table(bkG1, sep = "\t", file = "bkG1.txt")
    cat("The screened pairs are now in your working directory. ")
    list(bkG1 = bkG1, bkG2 = bkG2)
}

