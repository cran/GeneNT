getBM <- function(pG2, kG2)
{
   ##combine two screened pairs
   ppair <- pG2[,3:5]
   kpair <- kG2[,3:5]
   pair <- rbind(ppair, kpair)

   #find unsorted unique list of probsets
   index.us <- unique(rbind(as.matrix(pair[,1]), as.matrix(pair[,2])))
   #construct Boolean connectivity matrix to export to Pajek
   BM <- matrix(0, nrow(index.us), nrow(index.us)) 
   for(i in 1:nrow(pair))
   {
 	x <- which(index.us == as.character(pair[i,1])) 
 	y <- which(index.us == as.character(pair[i,2])) 
 	BM[x,y] <- 1
 	BM[y,x] <- 1
   }
   diag(BM) <- 1
   row.names(BM) <- index.us
   colnames(BM) <- index.us
   write.table(BM, sep = "\t", file = "BM.tsv")
   
   #the code belwo is obtained from http://vlado.fmf.uni-lj.si/pub/networks/pajek/howto/HowToR.htm
   #to save ordinary matrix in R to pajek compatible
   savematrix <- function(n,direct,twomode=1){
    if ((dim(n)[1] == dim(n)[2]) & (twomode!=2))
      { write(paste("*Vertices",dim(n)[1]), file = direct);
            write(paste(seq(1,length=dim(n)[1]),' "',rownames(n),
                  '"',sep=""), file = direct,append=TRUE);
            write("*Matrix", file = direct,append=TRUE);
            write(t(n),file = direct,ncolumns=dim(n)[1],
                  append=TRUE) }
    else
      { write(paste("*Vertices",sum(dim(n)),dim(n)[1]),
              file = direct);
            write(paste(1:dim(n)[1],' "',rownames(n),'"',sep=""),
                  file = direct,append=TRUE);
            write(paste(seq(dim(n)[1]+1,length=dim(n)[2]),' "',
                  colnames(n),'"',sep=""), file = direct,append=TRUE);
            write("*Matrix", file = direct, append=TRUE);
            write(t(n),file = direct, ncolumns=dim(n)[2],append=TRUE)}
      }   
   savematrix(BM,"BMPajek.mat")
}

