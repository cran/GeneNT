spclust <- function(p, pG2, kG2)
{
   ##combine two screened pairs
   ppair <- pG2[,3:5]
   kpair <- kG2[,3:5]
   pair <- rbind(ppair, kpair)

   #find unsorted unique list of probsets
   index.us <- unique(rbind(as.matrix(pair[,1]), as.matrix(pair[,2])))

   #calculate the connectivity list
   mdegree <- matrix(0, nrow(index.us), nrow(index.us)) 
   for(i in 1:nrow(pair))
   {
 	x <- which(index.us == as.character(pair[i,1])) 
 	y <- which(index.us == as.character(pair[i,2])) 
 	mdegree[x,y] <- 1
 	mdegree[y,x] <- 1
   }

   #order probsets according to connecitvity
   ll <- apply(mdegree, 1, sum)
   ll <- as.matrix(ll)
   degreedata <- cbind(index.us, ll)
   idx <- order(as.numeric(degreedata[,2]), decreasing = T)
   index <- degreedata[idx,1]
   index <- data.frame(index)

   M <- matrix(NA, nrow(index), nrow(index)) 
   for(i in 1:nrow(pair))
   {
 	x <- which(index == as.character(pair[i,1])) 
 	y <- which(index == as.character(pair[i,2])) 
        if(is.na(M[x,y])) 
  	{
   	   M[x,y] <- (1- abs(pair[i,3]))^p
           M[y,x] <- (1- abs(pair[i,3]))^p
 	} 
 	if(!is.na(M[x,y])) #choose the largest correlation, shortest path
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
   #mdegree is the distance matrix filled with SPs (SP distance matrix)
   Mdist <- z1$length
   colnames(Mdist) <- rcname
   row.names(Mdist) <- rcname

   #find the GSC.
   if(!any(is.na(Mdist)))
   { gsc <- Mdist }
   if(any(is.na(Mdist)))
   { 
      i = 1
      while(1)
      {
 	 if(any(is.na(Mdist[1:i, 1:i])))
 	 break
 	 i <- i + 1 
      }
      gsc <- Mdist[1:(i-1), 1:(i-1)]
   }
   
   if(any(is.na(gsc))) #check whether NA is mising  
   cat("GCC search error!")
   
   write.table(gsc, sep = "\t", file = "gsc.tsv")
   d <- as.dist(gsc)
   obj <- hclust(d, method = "single")
   write.table(obj$label[obj$order], sep = "\t", file = "label.tsv")   
   cat("The Shortest-Path distance matrix and labels have been written to the default directory!")   
   plot(obj, label = F, hang = 0, main = "Clustering with posterior (shortest-path) distance matrix", sub = "", xlab = "" )
}

