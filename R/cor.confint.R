cor.confint <- function (cor, N, alpha) 
{
        z <- atanh(cor)
        se <- 1/sqrt(N - 3)
        conf.int1 <- tanh(z-qnorm(1- alpha/2)*se)        
        conf.int2 <- tanh(z+qnorm(1- alpha/2)*se)        
	list(conf.int1 = conf.int1, conf.int2 = conf.int2)
}

