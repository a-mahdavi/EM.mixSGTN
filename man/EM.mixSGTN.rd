\name{EM.mixSGTN}
\alias{EM.mixSGTN}
\title{EM.mixSGTN}
\usage{
EM.mixSGTN(g=1, w=1, xi, s, b, la, nu=1, iter.max=200, tol=10^-6, get.init = TRUE, family="SGTN",group=F, calc.im=FALSE)
}
\description{
 Fit the mixures of SGTN distributions using EM-algorithm
 g: number of components to fit
 w: the vector of components probabilites.
 xi: the vector of location parameters.
 s: the vector of scale parameters.
 la: the vector of shape parameters.
 nu: the flatness parameters.
 family: distribution family to be used in fitting ("SGTN", "SGN").
 get.init: if TRUE, the initial values are generated.
 iter.max: the maximum number of iterations of the EM algorithm. Default= 100.
 tol: the covergence maximum error. Default= 10^-6.
 group: if TRUE it returns the group id of each observation
 calc.im: if TRUE, the information matrix is calculated and the standard errors are reported
 }
\examples{

 #  Example 1:
 # Simulating 100 samples from one component SGTN distribution:
	y <- r.mixSGTN(n=100, xi=5, s=2, la=3, b=3, nu=5)
 # EM output with specific initial values: 
 EM.mixSGTN(y, xi=4, s=2.5, la=1, b=4, nu=3, get.init=FALSE)
 # EM output without specific initial values: 
 EM.mixSGTN(y, get.init=TRUE)

 # Example 2:
 # Simulating 1000 samples from mixtures of SGTN distributions:
	y <-  r.mixSGTN(n=1000, w=c(.3,.7), xi=c(0,5), s=c(2,5), b=c(1,3), la=c(-3,3) , nu=4)
  # EM output without specific initial values: 
 EM.mixSGTN(y, g=2, get.init=TRUE)
   }

