
EM.mixSGTN <- function(y, g=1, w=1, xi, s, b, la, nu=2, iter.max=100, tol=10^-6, get.init = TRUE, family="SGTN",group=F, calc.im=FALSE){  
 # g: number of components to fit
 # w: the vector of components probabilites.
 # xi: the vector of location parameters.
 # s: the vector of scale parameters.
 # la: the vector of shape parameters.
 # nu: the flatness parameters.
 # family: distribution family to be used in fitting ("SGTN", "SGN").
 # get.init: if TRUE, the initial values are generated.
 # iter.max: the maximum number of iterations of the EM algorithm. Default= 100.
 # tol: the covergence maximum error. Default= 10^-6.
 # group: if TRUE it returns the group id of each observation
 # calc.im: if TRUE, the information matrix is calculated and the standard errors are reported
	begin <- proc.time()[3] 
	dSGT <- function(y, xi, s, b, la, nu){
		d <- b/(nu^(1/b)*beta(nu/b,1/b)*s*(1+abs(y-xi)^b/(nu*s^b))^((nu+1)/b))*pnorm(la*(y-xi)/s)
		d[which(d==0)]=10^-100
		return(d) }
	dmixSGT <- function(y, w, xi, s, b, la, nu){
			d <- 0 ; g <- length(w)
			for ( j in 1:g)
			d <- d + w[j]*dSGT(y, xi[j], s[j], b[j], la[j], nu)
			return(d) }
	skewness <- function (x, na.rm = FALSE) 
		{
    	if (is.matrix(x)) 
        apply(x, 2, skewness, na.rm = na.rm)
    	else if (is.vector(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        n <- length(x)
        (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  	  	}
    	else if (is.data.frame(x)) 
        sapply(x, skewness, na.rm = na.rm)
    	else skewness(as.vector(x), na.rm = na.rm)
		}
	n <- length(y)   ;     dif <- 1
        count <- 0 
	if (get.init == TRUE) {
               init <- kmeans(y, g,  algorithm="Hartigan-Wong")
               w <- init$size/n ; la <- NULL
               xi <- as.vector(init$centers)
               s <- sqrt(init$withinss/init$size) 
			for (k in 1:g)	
		   la[k] <- skewness(y[init$cluster==k])				
			b <- runif(g, 1,5)
		}
	 LL <- 1 
	while ((dif > tol) && (count <= iter.max)) {
	z.hat <- gam.hat <- matrix(0,n,g)
# E step
	for (j in 1:g){
	z.hat[,j] <- w[j]*dSGT(y,xi[j],s[j],b[j],la[j], nu)/dmixSGT(y, w, xi, s, b ,la, nu)
	y.xi <- y-xi[j] ; del<-la[j]/s[j]; del.y.xi <- del*y.xi; eta=abs(y.xi)/s[j]
	del.y.xi[which(del.y.xi< (-30))]=-30
	s1 <- (del.y.xi + dnorm(del.y.xi)/pnorm(del.y.xi))
	s2 <- ((1+nu)/(nu+eta^b[j]))
# MCE steps
	w[j] <- sum(z.hat[,j])/n
	y.xi[which(y.xi==0)]=10^-6
	xi[j] <- (sum(z.hat[,j]*(s2/(s[j]^b[j])*abs(y.xi)^(b[j]-2)+del^2)*y)-del*sum(z.hat[,j]*s1))/sum(z.hat[,j]*(s2/(s[j]^b[j])*abs(y.xi)^(b[j]-2)+del^2))
	y.xi <- y-xi[j]
	s[j] <- (sum(z.hat[,j]*s2*abs(y.xi)^b[j])/(sum(z.hat[,j])))^(1/b[j])
	eta=abs(y.xi)/s[j]
	del <- sum(z.hat[,j]*s1*y.xi)/sum(z.hat[,j]*y.xi^2)
	la[j]=del*s[j]
	if( family=="SGTN"){
	cml <- optim(c(b[j],nu),function(x){
		b <- x[1] ; nu <- x[2]
		-sum(z.hat[,j]*log(dSGT(y,xi[j],s[j],b,la[j],nu)))
		},method="L-BFGS-B",lower=c(.1,.1),upper=c(30,100))$par
	b[j] <- cml[1] ; nu <- cml[2] }
	else {
	nu <- 100
	b[j] <- optim(b[j],function(x){
		b <- x[1] 
		-sum(z.hat[,j]*log(dSGT(y,xi[j],s[j],b,la[j],nu)))
		},method="L-BFGS-B",lower=c(.1,.1),upper=c(30,100))$par 
	} }
	LL.new <- sum(log(dmixSGT(y,w,xi,s,b,la,nu))) # log-likelihood function
	count <- count +1 
	dif <- abs(LL.new/LL-1)
	}
	if( calc.im ){
	if( family=="SGTN"){
	IM <- 0
	g <- length(w); sw <-sxi<-ss<-sla<-sb<-snu<- NULL ; nu=rep(nu,g)
	for( i in 1:length(y)){
	for(j in 1:(g-1))
	sw[j] <- z.hat[i,j]/w[j]-z.hat[i,g]/w[g]
	for(j in 1:g){
	z.hat[i,j] <- w[j]*dSGT(y[i],xi[j],s[j],b[j],la[j], nu[1])/dmixSGT(y[i], w, xi, s, b ,la, nu[1])
	s1 <- (la[j]*(y[i]-xi[j])/s[j] + dnorm(la[j]*(y[i]-xi[j])/s[j])/pnorm(la[j]*(y[i]-xi[j])/s[j]))
	if( is.na(s1) || s1==Inf) 
	s1 <- (la[j]*(y[i]-xi[j])/s[j] + dnorm(la[j]*(y[i]-xi[j])/s[j])/pnorm(-30))
	s2 <- ((1+nu[j])/(nu[j]+(abs(y[i]-xi[j])/s[j])^b[j]))
	kappa <- digamma((nu[j]+1)/b[j])+log((nu[j]+(abs(y[i]-xi[j])/s[j])^b[j])/b[j])
	zeta <- abs(y[i]-xi[j])/s[j]
	sxi[j] <- z.hat[i,j]/s[j]*(s2*zeta^(b[j]-1)*sign(y[i]-xi[j])-la[j]*(s1-la[j]*zeta))
	ss[j] <- z.hat[i,j]/s[j]*(s2*zeta^b[j]+la[j]^2*zeta^2-la[j]*s1*(y[i]-xi[j])/s[j]-1)
	sla[j] <- z.hat[i,j]/s[j]*(y[i]-xi[j])*(s1-la[j]*(y[i]-xi[j])/s[j])
	sb[j] <- z.hat[i,j]/b[j]^2*(digamma(1/b[j])+nu[j]*digamma(nu[j]/b[j])+log(b[j])+b[j]
         +s2*zeta^b[j]*(1-b[j]*log(abs(y[i]-xi[j])))+nu[j]*s2+kappa-nu[j]*log(nu[j]/b[j])
	   -nu[j]/b[j]^2-1 )
	snu[j] <- z.hat[i,j]/b[j]*(log(nu[j]/b[j])+1-digamma(nu[j]/b[j])-s2)
	}
	IM <- IM + c(sw,sxi,ss,sb,sla,sum(snu))%*%t(c(sw,sxi,ss,sb,sla,sum(snu)))
	}
	se <- sqrt(diag(solve(IM)))
	}
	else{
	IM <- 0
	g <- length(w); sw <-sxi<-ss<-sla<-sb<-snu<- NULL ; nu=rep(100,g)
	for( i in 1:length(y)){
	for(j in 1:(g-1))
	sw[j] <- z.hat[i,j]/w[j]-z.hat[i,g]/w[g]
	for(j in 1:g){
	z.hat[i,j] <- w[j]*dSGT(y[i],xi[j],s[j],b[j],la[j], nu[1])/dmixSGT(y[i], w, xi, s, b ,la, nu[1])
	s1 <- (la[j]*(y[i]-xi[j])/s[j] + dnorm(la[j]*(y[i]-xi[j])/s[j])/pnorm(la[j]*(y[i]-xi[j])/s[j]))
	if( is.na(s1) || s1==Inf) 
	s1 <- (la[j]*(y[i]-xi[j])/s[j] + dnorm(la[j]*(y[i]-xi[j])/s[j])/pnorm(-30))
	s2 <- ((1+nu[j])/(nu[j]+(abs(y[i]-xi[j])/s[j])^b[j]))
	kappa <- digamma((nu[j]+1)/b[j])+log((nu[j]+(abs(y[i]-xi[j])/s[j])^b[j])/b[j])
	zeta <- abs(y[i]-xi[j])/s[j]
	sxi[j] <- z.hat[i,j]/s[j]*(s2*zeta^(b[j]-1)*sign(y[i]-xi[j])-la[j]*(s1-la[j]*zeta))
	ss[j] <- z.hat[i,j]/s[j]*(s2*zeta^b[j]+la[j]^2*zeta^2-la[j]*s1*(y[i]-xi[j])/s[j]-1)
	sla[j] <- z.hat[i,j]/s[j]*(y[i]-xi[j])*(s1-la[j]*(y[i]-xi[j])/s[j])
	sb[j] <- z.hat[i,j]/b[j]^2*(digamma(1/b[j])+nu[j]*digamma(nu[j]/b[j])+log(b[j])+b[j]
         +s2*zeta^b[j]*(1-b[j]*log(abs(y[i]-xi[j])))+nu[j]*s2+kappa-nu[j]*log(nu[j]/b[j])
	   -nu[j]/b[j]^2-1 )
	}
	IM <- IM + c(sw,sxi,ss,sb,sla)%*%t(c(sw,sxi,ss,sb,sla))
	}
	se <- sqrt(diag(solve(IM)))
	}
	}
	if( family=="SGTN"){
	aic <- -2 * LL.new + 2 * (4*g+g)
	bic <- -2 * LL.new + log(n) * (4*g+g)
	edc <- -2 * LL.new + 0.2*sqrt(n) * (4*g+g)
	}
	else {
	aic <- -2 * LL.new + 2 * (4*g+g-1)
	bic <- -2 * LL.new + log(n) * (4*g+g-1)
	edc <- -2 * LL.new + 0.2*sqrt(n) *(4*g+g-1)
	}
	end <- proc.time()[3]
	time <- end-begin
	obj.out <- list(family=family,w=w, xi=xi, sigma=s, b=b, la=la , nu=nu[1], loglik=LL.new, aic=aic, bic=bic, edc=edc, iter=count,elapsed=as.numeric(time))
	if (group)
	obj.out$group <- apply(z.hat, 1, which.max)
	if (calc.im)
	obj.out$std <- se
	obj.out
	}


