BEST.pearson <- function(gene.name, Q, cormin, method = c("simple", "jackknife"))
{
    if (class(dat) != "matrix")
    dat <- as.matrix(dat)
    qy <- which(names(dat[,0][1:nrow(dat)])==gene.name)
    ##########
    method = match.arg(method)
    p.list <- array(NA, nrow(dat))
    cor.list <- array(NA, nrow(dat))
    if (method == "simple")
    {
        for (j in 1:nrow(dat))
        {
        g<- cor.test(dat[qy,], dat[j,], alternative = "two.sided", method = "pearson")
        cor.list[j] <- cor(dat[qy,], dat[j,])
        p.list[j] <- g$p.value
        }
    }
    if (method == "jackknife")
    {
        cor.temp <- array(NA, ncol(dat))
        for (j in 1:nrow(dat))
        {
            for (k in 1:ncol(dat))
            {
            if((sd(dat[qy,-k])*sd(dat[j,-k]))==0)
            cor.temp[k] <- 0
            if((sd(dat[qy,-k])*sd(dat[j,-k]))!=0)
            cor.temp[k] <- cor(dat[qy,-k], dat[j, -k], method = "pearson")
            }
        cor.list[j] <- median(cor.temp) 
        g<- cor.test(dat[qy,], dat[j,], alternative = "two.sided", method = "pearson")
        p.list[j] <- g$p.value
        }
    }
    ####calculate and rank p-values for all pairs
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

    ###calculate q values (FDR p-values)
    fdr.out <- fdr.control(p.list, Q)
    q.list <- fdr.out$qvalues
    q.results <- cbind(p.results, q.list)

    ###use stepdown procedure to generate G1 set controlling FDR
    G1.results <- q.results[(q.results[,7] < Q),]
    sort.idx <- order(-abs(G1.results[, 5]))
    G1.results <- G1.results[sort.idx, ]
    G1 <- nrow(G1.results)
    pq <- G1/length(p.list)

    ###construct confidence intervals for G1
    lower <- array(NA, G1)
    higher <- array(NA, G1)
    if(G1 == 0)
    warning("Stage I screening returns no results. The Q value may be set to conservative (low)!")		
    for(i in 1:G1)
    {
        x <- G1.results[i,1]
        y <- G1.results[i,2]
        x.row <- dat[x,]
        y.row <- dat[y,]
        alpha = Q*pq
        g <- cor.test(x.row, y.row, alternative = "two.sided", method = "pearson", conf.level = 1- alpha/2)
        lower[i] <- g$conf.int[1]
        higher[i] <- g$conf.int[2]
     }
    lower <- as.matrix(lower)
    lower <- data.frame(lower)
    higher <- as.matrix(higher)
    higher <- data.frame(higher)
    bpG1 <- cbind(G1.results, lower, higher) 
    indxg2 <- apply(bpG1,1,function(x) ifelse(((as.numeric(x[8]) > cormin) | (as.numeric(x[9]) < -cormin)), indxg2<- T, indxg2 <- F))
    bpG2 <- bpG1[indxg2,]
    write.table(bpG2, sep = "\t", file = "bpG2.txt")
    write.table(bpG1, sep = "\t", file = "bpG1.txt")
    cat("The screened pairs are now in your working directory. ")
    list(bpG1 = bpG1, bpG2 = bpG2)

}

