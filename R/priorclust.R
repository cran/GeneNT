priorclust <- function(p)
{
   #dat <- read.table(file.name, h = T,row.names = 1) 
   C <- cor(t(dat))
   diag(C) <- 0
   D <- (1-abs(C))^p
   diag(D) <- 0
   write.table(D, sep = "\t", file = "D.tsv")
   d <- as.dist(D)
   obj <- hclust(d)
   plot(obj, label = F, hang = 0, main = "Clustering with prior distance matrix", xlab = "" )
   y <- identify(obj, MAXCLUSTER = 50)
}

