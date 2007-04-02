kendall.confint <- function (x, y, alpha) 
{
	##the followed code is adopted in large part from http://www.stat.umn.edu/geyer/5601/examp/corr.html
        signs <- sign(outer(x, x, "-") * outer(y, y, "-"))
	tau <- mean(signs[lower.tri(signs)])
	#tau <- cor(x,y,method = "kendall")
	cvec <- apply(signs, 1, sum)
	n <- length(cvec)
	sigsq <- (2 / (n * (n - 1))) * (((2 * (n - 2)) / (n * (n - 1))) * var(cvec) + 1 - tau^2)
	zcrit <- qnorm(1 - alpha/2)
	conf.int <- tau + c(-1, 1) * zcrit * sqrt(sigsq)
	conf.int
}

