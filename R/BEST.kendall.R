BEST.kendall <- function(gene.name, Q, cormin)
{
    if (class(dat) != "matrix")
    dat <- as.matrix(dat)    
    qy <- which(names(dat[,0][1:nrow(dat)])==gene.name)
    p.list <- array(NA, nrow(dat))
    cor.list <- array(NA, nrow(dat))
    for (j in 1:nrow(dat))
    {
        g<- cor.test(dat[qy,], dat[j,], alternative = "two.sided", method = "kendall", exact = F)
        cor.list[j] <- cor(dat[qy,], dat[j,], method = "kendall")
        p.list[j] <- g$p.value
    }
    ###Calculate p-values for all pairs
    p.name1 <- rep(names(dat[,0][qy]),dim(dat)[1])
    p.name2 <- names(dat[,0][1:dim(dat)[1]])
    p.name <- cbind(p.name1, p.name2)
    colnames(p.name) <- c("gene1", "gene2")
    p.name <- data.frame(p.name)
    p.index1 <- rep(qy, dim(dat)[1])  
    p.index2 <- seq(1:dim(dat)[1])
    p.index <- cbind(p.index1, p.index2)
    colnames(p.index) <- c("index1", "index2")
    p.results <- cbind(p.index, p.name, cor.list, p.list)
    ###calculate q value
    fdr.out <- fdr.control(p.list, Q)
    q.list <- fdr.out$qvalues
    q.results <- cbind(p.results, q.list)
    ###use stepdown procedure to generate G1 set controlling FDR
    G1.results <- q.results[(q.results[,7] < Q),]
    sort.idx <- order(-abs(G1.results[, 5]))
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
        x <- G1.results[i,1] #get indices
        y <- G1.results[i,2]
        x <- dat[x,]
        y <- dat[y,]
        g <- kendall.confint(x,y,alpha)
        lower[i] <- g[1]
        if(g[2] >1)
            higher[i] <- 1
        if(g[2] <=1)
            higher[i] <- g[2]
    }
    lower <- as.matrix(lower)
    lower <- data.frame(lower)
    higher <- as.matrix(higher)
    higher <- data.frame(higher)    
    bkG1 <- cbind(G1.results, lower, higher) 
    indxg2 <- apply(bkG1,1,function(x) ifelse(((as.numeric(x[8]) > cormin) | (as.numeric(x[9]) < -cormin)), indxg2<- T, indxg2 <- F))
    bkG2 <- bkG1[indxg2,]
    write.table(bkG2, sep = "\t", file = "bkG2.txt")
    write.table(bkG1, sep = "\t", file = "bkG1.txt")
    cat("The screened pairs are now in your working directory. ")
    list(bkG1 = bkG1, bkG2 = bkG2)
}

