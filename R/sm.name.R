sm.name <- function (m) 
{
    l <- nrow(m)
    index1 <- rep(NA, l * (l - 1)/2)
    index2 <- rep(NA, l * (l - 1)/2)
    k <- 1
    for (i in 1:(l - 1)){ 
	for (j in (i + 1):l) {
        index1[k] <- names(m[i])
        index2[k] <- names(m[j])
        k <- k + 1
    }
    }
    return(cbind(index1, index2))
}

