tdclust <- function(p)
{
   C <- cor(t(dat))
   diag(C) <- 0
   D <- (1-abs(C))^p
   diag(D) <- 0
   write.table(D, sep = "\t", file = "D.tsv")
   d <- as.dist(D)
   obj <- hclust(d)
   plot(obj, label = FALSE, hang = 0, main = "Traditional clustering", sub = "", xlab = "" )
}

