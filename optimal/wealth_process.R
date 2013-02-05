### Properties of the alpha investing processes

	
plot.alpha <- function(alpha,range=c(0,0.01)) {
	w <- seq(0,1,0.01);
	plot(w,sapply(w,function(x){alpha(x,1)}),col="red", type="l", ylim=range);
	lines(w,sapply(w,function(x){alpha(x,2)}),col="green");
	lines(w,sapply(w,function(x){alpha(x,3)}),col="blue");
	}

######################################################################
###   
###   omega is global constant set by these functions
###
      omega <- function() { 0.5 }
###
######################################################################

K.scale <- 4;

###  Harmonic  ###
#
#     Need the constant larger than a multiple of 1/p.reject or wealth 
#     will go negative.  This produces a large value for the position x.
#     With a very large x, then you spend next to nothing until the process
#     finds the right level (i.e., would have lower power).
#

alpha.har <- function(w, scale=K.scale) {
	K.har <- log(1000);
	Finv  <-  function(F) { exp(K.har-F) };
	f     <-  function(x) { 1/x };
	scale*f(Finv(w/scale))
}
#	F     <-  function(x) { K.har - log(x) };


###  Universal ### (note you can solve for x < 1 to handle high wealth)

alpha.univ <- function(w, scale=K.scale) {
	K.univ <- 1/c(2.10974, 1.06906, 0.79288)[d.univ <- 1];
	Finv <- function(F) { exp(K.univ/F)-d.univ };
	f    <- function(x) { K.univ/((x+d.univ)*log(x+d.univ)^2) };
	scale*f(Finv(w/scale))
}
#	F <- function(x) { (K.univ/log(x+d.univ)) },


###  Geometric ### (note you can solve for x < 1 to handle high wealth)

alpha.geo <- function(w, scale=K.scale){
	psi <- 0.05;
	psi * w
}
#	Finv <-function(F) { z <- log(1-psi); log(F*(psi-1)*z/psi)/z };
#	f <-	function(x) { psi*(1-psi)^(x-1) };
#	scale*f(Finv(w/scale))
#	F <- 	function(x) { psi*(1-psi)^(x-1)/log(1-psi) },


#####################################################################
# All had a problem that |f| > F for some x < 1.  That means they
# will try to spend more than they have, which is not good.  
#   Remedy: bound derivative by omega.
#
#####################################################################

check <- function() {
	# --- check inverse functions
	Finv(F(10))
	F(Finv(0.3))
	# --- plot to see wealth function F and bids (neg ok for geometric)
	plot.alpha(alpha.univ,c(0,0.3))

	F(x)                  # wealth at x
	Finv(F(x)) - x        # should be 0

	f(x)                  # amount spent
	F(x) + f(x)           # wealth after spending
	
	Finv(F(x) + f(x))     # next position if do not reject
	Finv(F(x) + 0.5)      # new position if reject

	rm(Finv, F, f, parms)
	}
	
#####################################################################

simulate <- function(nSteps) {
	# --- conditions
	p.reject <<- 0.05; 
	p.values <- rbinom(nSteps,1,1-p.reject) # 0/1 pvalues

	wealth <- matrix(0, nrow=nSteps, ncol=2);
	colnames(wealth) <- c("geo","univ");

	wealth[1,] <- rep(0.5	,2)
	for(i in 2:nSteps)	{ 
		a <- alpha.geo(w <-wealth[i-1,"geo"])
		wealth[i,"geo"] <- ifelse(p.values[i]<a, w+omega(), w-a)
		a <- alpha.univ(w <-wealth[i-1,"univ"])
		wealth[i,"univ"] <- ifelse(p.values[i]<a, w+omega(), w-a)
		}
	wealth
}


# --- where does it stay
draw.hist <- function() {
	hist(log(wealth[,"geo"]), 
		breaks=100, main="Geometric", xlab="Wealth",
		sub=paste("p=",p.reject,"  omega=",omega())
		)
	hist(log(wealth[,"univ"]), 
		breaks=100, main="Universal", xlab="Wealth",
		sub=paste("p=",p.reject,"  omega=",omega())
		)
	}

see.results <- function() {

	nSteps <- 2000; 	i <- 1:nSteps

	wealth <- simulate(nSteps)

	idx <- "geo"; plot(i, wealth[i,idx], main=paste("Wealth(",idx,")"), type="p")
	idx <- "univ"; points(i, wealth[i,idx], pc=23)
	 
	idx <- "geo"; i <- 1:500
	plot(  i[-1], -diff(wealth[i,idx]), main="Changes", type="l", ylim=c(0,0.05))
	idx <- "univ"; lines(i[-1], -diff(wealth[i,idx]), col="red")

	par(mfrow=c(3,1))
		draw.hist(simres[,1] , parms.geo );
		draw.hist(simres[,2] , parms.univ );
		draw.hist(simres[,3] , parms.har);
	par(mfrow=c(1,1))

	with(parms.univ, F(0.54462))

# ---- harmonic, constant K set to 1000, x0=606.5
# 0.100, ½     21.6 -0.0462
# 0.050        45.6 -0.0219
# 0.025        96.8 -0.0103
# 0.010       272.4 -0.00367

# ---- harmonic, K = 10000, x0 = 6065, omega = ½ 
# 0.1000       21.5  -0.047     -0.055
# 0.0500       47.4  -0.021     -0.025
# 0.0100      251    -0.0042    -0.0048
# 0.0050      489    -0.0021    -0.0023
# 0.0025      963    -0.0012    -0.0012

# ----  harmonic, K = 10000,  omega = ¼ 
# 0.1000       37.3  -0.027     -0.029
# 0.0500       84.8  -0.012     -0.013  
# 0.0100      423    -0.0024    -0.0024

# ----  Geometric, psi = 0.1, omega = ½ 
# 0.10          8.3  -0.046     -0.054
# 0.05         16.7  -0.019     -0.028
# 0.01         77.9  -0.0003    -0.005

# --- sequence of states
	plot(x, pch=19, cex=0.25)

	plot(x, xlim=c(4000,6000), ylim=c(0,400), pch=19, cex=0.25)
	

# --- what's the process look like, superimposed on wealth function
	plot(F,0, 1.1*max(x))
	points(x,sapply(x,F), col="blue")

	plot(F,0, 250); points(x,sapply(x,F), col="blue", pch=19, cex=0.25)

# --- comparisons … similar durations since same arrival process
	par(mfrow=c(2,1))
	  plot(x.geo, pch=19, cex=0.25)
	  plot(x.har, pch=19, cex=0.25)
	par(mfrow=c(1,1))

	par(mfrow=c(2,1))
	  bid <- -f.geo(x.geo); plot(bid[(0<bid)&(bid<0.05)], pch=19, cex=0.25, main="Geometric")
 	 bid <- -f.har(x.har); plot(bid[(0<bid)&(bid<0.05)], pch=19, cex=0.25, main="Harmonic")
	par(mfrow=c(1,1))

	par(mfrow=c(2,1))
	  bid <- -f.geo(x.geo); hist(bid[(0<bid)&(bid<0.025)],breaks=100,  main="Geometric",
	  		sub=paste("SD = ",round(sd(bid),4)))
 	 bid <- -f.har(x.har); hist(bid[(0<bid)&(bid<0.025)],breaks=100,  main="Harmonic",
  			sub=paste("SD = ",round(sd(bid),4)))
	par(mfrow=c(1,1))

	par(mfrow=c(2,1))
 		bid <- -f.geo(x.geo); hist(log10(bid[(0<bid)&(bid<1)]), breaks=100, main="Geometric")
  		bid <- -f.har(x.har); hist(log10(bid[(0<bid)&(bid<1)]), breaks=100, main="Harmonic")
	par(mfrow=c(1,1))
}




