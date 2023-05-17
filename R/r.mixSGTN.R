r.mixSGTN <- function(n , w=1, xi, s, b, la, nu){
	rgnorm <- function (n, mu = 0, alpha = 1, beta = 1) 
	{
    if (alpha <= 0 | beta <= 0) {
        cat("Not defined for negative values of alpha and/or beta.\n")
        return(rep(NaN, n))
    }
    lambda <- (1/alpha)^beta
    unifs <- runif(n)
    scales <- qgamma(unifs, shape = 1/beta, scale = 1/lambda)^(1/beta)
    return(scales * ((-1)^rbinom(n, 1, 0.5)) + mu)
	}
	r.SGTN <- function(n , xi, s, b, la, nu){
		y <- 0 
		for(i in 1:n){
		V <- rgamma(1,shape=nu/b,rate=nu/b)^(-1/b)*rgnorm(1, mu = 0, alpha = 1, beta = b)*b^(1/b)
		x1 <- rnorm(1)
		if( x1 < la*V )
		y[i] <- V
		else
		y[i] <- -V
		}
		return(xi + s*y)
		}
	g <- length(xi) ; y <-NULL ; Z <- rmultinom(n,size=1,prob=w); n <- rowSums(Z)
	for (j in 1:g)
	y <- c(y, r.SGTN(n[j], xi[j], s[j],  b[j], la[j], nu ) )
		return(y)
		}