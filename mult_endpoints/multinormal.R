##  Computing multivariate normal distribution

#########################
#  Note:
#    Need to check that 0,0,0... for mean gives right rejection region
#
#########################

library(mvtnorm)

#--------------------------------------------------------------------------------


# --- define correlations; these stay fixed for all of the tests

corr <- matrix(data=c(1,.25,.31,.24,  .25,1,.42,.43,   .31,.42,1,.67,   .24,.43,.67,1),nrow=4);

v3 <- matrix(data=c(1,-.25,-.31,  -.25,1,.42,   -.31,.42,1),nrow=3);

sig1 <- 1
sig2 <- corr[1:2,1:2]
sig3 <- corr[1:3,1:3]
sig4 <- corr


#--------------------------------------------------------------------------------
#  Compare integrals to those in MMa

tau1 <- 1.96
tau2 <- 2.49

g3 <- function(tau3) {
 pmvnorm(c(0,0,0),sig3,lower=c(tau1,-Inf,tau3),upper=c(Inf,tau2,Inf))/
   pmvnorm(c(0,0),sig2,lower=c(tau1,-Inf),     upper=c(Inf,tau2))
   }
   
# Why on earth are these so different?  And why does it change each time its called?
round(g3(2), 6)
round( mapply(g3,0:4), 6)

# --- do as in MMa by flipping signs so its ratio of CDFs

g3s <- function(tau3) {
 pmvnorm(c(0,0,0),v3,         lower=c(-Inf,-Inf,tau3),  upper=c(-tau1,tau2,Inf))/
   pmvnorm(c(0,0),v3[1:2,1:2],lower=c(-Inf,-Inf),       upper=c(-tau1,tau2)  )
   }
g3s(2)
round( mapply(g3s,0:4), 6)

# --- compute the integral with a bit of averaging replications

rasnorm <- function(m,v,lower,upper) {
	results <- rep(NA,50);
	for (i in 1:length(results)) { results[i] <- pmvnorm(m,v,lower=lower,upper=upper); }
    mean(results)
    }
    
g3r <- function(tau3) {
 rasnorm(c(0,0,0),v3,         lower=c(-Inf,-Inf,tau3),  upper=c(-tau1,tau2,Inf))/
   rasnorm(c(0,0),v3[1:2,1:2],lower=c(-Inf,-Inf),       upper=c(-tau1,tau2)  )
   }
g3r(2)
       

#--------------------------------------------------------------------------------
# initial scalar case

tau1 <- qnorm(1-0.025); tau1                  # TAU1

1-pnorm(tau1)


#--------------------------------------------------------------------------------
# 2 terms
#                H1 rejects -->  mu1 ≥ 0

# --- first look at placement of the means
mu1 <- seq(  0, 1, .01)
mu2 <- seq(-.6, 0, .05)

# set some initial location for tau2 from MMa 
tau2 <- 2.49

# --- plot probabilities of rejection region for this tau over field of choices for mu1 and mu2
#     looks like max is at the origin
p <- matrix(0,nrow=length(mu1),ncol=length(mu2))
rownames(p) <- paste("mu1=",mu1); colnames(p) <- paste("mu2=",mu2)
for(i in 1:length(mu1)) for(j in 1:length(mu2))
  p[i,j]<-pmvnorm((c(mu1[i],mu2[j])),sig2,lower=c(tau1,tau2),upper=c(Inf,Inf))/(1-pnorm(tau1,mean=mu1[i]))

round(p,4)
contour(mu1,mu2,p)

# --- plot as function of tau2 at origin
f <- function(t) { pmvnorm(c(0,0),sig2,lower=c(tau1,t),upper=c(Inf,Inf))/(1-pnorm(tau1,mean=0)) }
x <- seq(2.7,2.9,0.01)
plot(x,mapply(f,x),type="b"); abline(h=0.0125,col="red")


# --- results for 2-dim, using my modified version of pmvnorm
mu2 <- c(0,0)
tau2 <- 2.49
rasnorm(mu2,sig2,lower=c(tau1,tau2),upper=c(Inf,Inf))/(1-pnorm(tau1,mean=mu2[1]))




#--------------------------------------------------------------------------------
# 3 terms
#                   H1 rejects, H2 does not
#
#   Area is different from MMa integration?  Use my averaged version, still very different

mu3 <- c(0,0,0)

tau3 <- 2.60

g3 <- function(tau){
 1-rasnorm(c(0,0,0),v3,lower=c(-Inf,-Inf,-Inf),   upper=c(-tau1,tau2,tau))/
   rasnorm(c(0,0),v3[1:2,1:2],lower=c(-Inf,-Inf), upper=c(-tau1,tau2))}

t <- seq(2.5,2.6,0.005)
plot(t,mapply(g3,t), type="b"); abline(h=0.025)

g3(2.557)

g3(2.6)

# --- see how it depends on the mean
g3m <- function(m,j, tau3=2.557) {
	mv <- c(0,0);
	mv[j] <- m;
	rasnorm(c(mv,0),sig3,lower=c(tau1,-Inf,tau3),upper=c(Inf,tau2,Inf))/
		rasnorm(mv,sig2,lower=c(tau1,-Inf),     upper=c(Inf,tau2)) }
g3m(0,1)

# it does indeed grow the the value of mu2
g3m(-1,2)



#--------------------------------------------------------------------------------
# 4 terms
#                  at point of numerical instability
mu4 <- c(0,0,0,0)

tau4 <- 1.3

g4 <- function(tau){
 pmvnorm(c(0,0,0,0),sig4,lower=c(tau1,-Inf,tau3,tau),upper=c(Inf,tau2,Inf,Inf))/
   pmvnorm(c(0,0,0),sig3,lower=c(tau1,-Inf,tau3),     upper=c(Inf,tau2,Inf))}

g4(3.1)

t <- seq(3.09,3.10,0.001)
plot(t,mapply(g4,t), type="b"); abline(h=0.075)

g4(3.0945)

# --- see how the probability depends on the mean
g4m <- function(m,j, tau4=3.0945) {
	mv <- c(0,0,0);
	mv[j] <- m;
	pmvnorm(c(mv,0),sig4,lower=c(tau1,-Inf,tau3,tau4),upper=c(Inf,tau2,Inf,Inf))/
		pmvnorm(mv,sig3,lower=c(tau1,-Inf,tau3),     upper=c(Inf,tau2,Inf)) }

mv <- seq(0,1,.1)

#     first mean is positive since rejected; plot has negative slope
plot( mv, mapply(function(m){g4m(m,1)},mv))

#     second mean is negative since did not reject; plot has positive slope
plot(-mv, mapply(function(m){g4m(m,2)},mv))

#     third mean is positive since rejected
#plot( mv, mapply(function(m){g4m(m,3)},mv))






# way tedious to find all of these, so rely on a numerical optimizer to explore space...

# function to minimize is squared deviation of probability from 0.9875; numerical error building up
m1.0 <- 4.76; m3.0 <- 8.7; tau4.0 <- 2.5
f <- function(par) { (0.9875-pmvnorm(c(par[1],0,par[2],0),sig4,lower=c(tau1,-Inf,tau3,-Inf),upper=c(Inf,tau2,Inf,par[3])))^2 }

opt <- optim(c(m1.0,m3.0,tau4.0),f); opt
f(opt$par)
