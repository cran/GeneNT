pcor.confint <- function (pcor, kappa, alpha) 
{
        z <- atanh(pcor)
        se <- 1/sqrt(kappa - 2)
        conf.int1 <- tanh(z-qnorm(1- alpha/2)*se)        
        conf.int2 <- tanh(z+qnorm(1- alpha/2)*se)        
	list(conf.int1 = conf.int1, conf.int2 = conf.int2)
}

