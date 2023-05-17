\name{r.mixSGTN}
\alias{r.mixSGTN}
\title{r.mixSGTN function}
\usage{
r.mixSGTN(n ,w, xi, s, la, nu=NULL, family="SMSTN.mix")
}
\description{
 Generating random samples from mixtures of SMSMN distributions 
 n: number of samples.
 w: the vector of probability of each component.
 xi: the list of vector of location parameters.
 s: the list of cov-variance matrices.
 la: the list of vector of shape parameters.
 nu: the list of vector of flatness parameters.
 family: distribution family to be used in fitting ("SMSTN.mix", "SMSSLN.mix", "SMSCN.mix").
	}
\examples{
 # Example 1:
 # Simulating 100 samples from one component SGTN distribution:
 	y <- r.mixSGTN(n=100, xi=5, s=2, la=3,	nu=5)

 # Example 2:
 # Simulating 1000 samples from mixtures of SGTN distributions:
	y <-  r.mixSGTN(n=1000, w=c(.3,.7), xi=c(0,5), s=c(2,5), b=c(1,3), la=c(-3,3) , nu=4)

   }


