ncclust <- function(p, pG2, kG2)
{
   ##combine two screened pairs
   ppair <- pG2[,3:5]
   kpair <- kG2[,3:5]
   pair <- rbind(ppair, kpair)
   index <- unique(rbind(as.matrix(pair[,1]), as.matrix(pair[,2])))
   M <- matrix(NA, nrow(index), nrow(index)) 
   
   g <- ftM2graphNEL(unique(as.matrix(pair[,1:2])), W=NULL, V=NULL, edgemode="undirected")
   CC <- connectedComp(g)
   gcc.list <- CC$"1"
   
   for(i in 1:nrow(pair))
   {
 	x <- which(index == as.character(pair[i,1])) 
 	y <- which(index == as.character(pair[i,2])) 
        if(is.na(M[x,y])) 
  	{
   	   M[x,y] <- (1- abs(pair[i,3]))^p
           M[y,x] <- (1- abs(pair[i,3]))^p
 	} 
 	#choose the largest correlation, shortest path
	if(!is.na(M[x,y]))  	
	{ 
   	   M[x,y] <- min((1- abs(pair[i,3]))^p, M[x,y])   
   	   M[y,x] <- M[x,y]
 	}     
   }
   rcname <- as.matrix(index)
   colnames(M) <- rcname
   row.names(M) <- rcname
   diag(M) <- 0
   
   z1 <- allShortestPaths(M)
   Mdist <- z1$length
   colnames(Mdist) <- rcname
   row.names(Mdist) <- rcname
   ##extract GCC
   gcc <- Mdist[gcc.list, gcc.list]

   d <- as.dist(gcc)
   obj <- hclust(d)
   write.table(obj$label[obj$order], sep = "\t", file = "label.tsv")   
   cat("The Network constrained distance matrix and dendrogram labels have been written to the default directory!")   
   plot(obj, label = F, hang = 0, main = "Network constrained clustering", sub = "", xlab = "" )
}

