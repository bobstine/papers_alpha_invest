### Properties of the alpha investing processes

# --- recursion...get next position from prior
#       p is rejection probability in Bayes model

build.recursion <- function(parms, p, omega) {
	function (x) { parms$Finv(parms$F(x) + ifelse(rbinom(1,1,p), omega, parms$f(x))) }
	}
	
plot.parms <- function(parms, range, w=.5) {
	x0 <- parms$Finv(w);
	plot(parms$F, range[1], max(x0+1,range[2]), main=parms$name, sub=paste("X0 =",round(x0,2)))
	lines(rep(x0,2), c(0,parms$F(x0)), col="red")
	}

######################################################################
###   
###   omega is global constant set by a function
###
      omega <- function() { 0.5 }
###
######################################################################

###  Harmonic  ###
#
#     Need the constant larger than a multiple of 1/p.reject or wealth 
#     will go negative.  This produces a large value for the position x.
#     With a very large x, then you spend next to nothing until the process
#     finds the right level (i.e., would have lower power).
#

K.har <- log(1000);

parms.har <- list(
	name = "Harmonic",
	F=      function(x) { K.har - log(x)}, 
	Finv =  function(F) { exp(K.har-F) },
	f =     function(x) { pmax(-omega(),-1/x) }
)


###  Universal ### (note you can solve for x < 1 to handle high wealth)
#
#   Choose the starting point to match up to the corresponding geometric
#

d.univ <- 3;
K.univ <- 1/c(2.10974, 1.06906, 0.79288)[d.univ]

parms.univ <- list(
	name = "Universal",
	F = 	function(x) { K.univ/log(x+d.univ) },
	Finv = function(F) { exp(K.univ/F)-d.univ },
	f = 	function(x) { pmax(-omega(),-K.univ/((x+d.univ)*log(x+d.univ)^2)) }
)


###  Geometric ### (note you can solve for x < 1 to handle high wealth)

psi <- 0.05
parms.geo <- list(
	name = "Geometric",
	F = 	function(x) { -psi*(1-psi)^(x-1)/log(1-psi) },
	Finv = function(F) { z <- log(1-psi); log(F*(psi-1)*z/psi)/z },
	f = 	function(x) { pmax(-omega(),-psi*(1-psi)^(x-1)) }
)


#####################################################################
# All have a problem that |f| > F for some x < 1.  That means they
# will try to spend more than they have, which is not good.  
#   Remedy: bound derivative by omega.
#
#####################################################################

check <- function() {
	# --- check inverse functions
	with( parms.har,  	Finv(F(10)))
	with( parms.har,  	F(Finv(0.3)))
	# --- plot to see wealth function F and bids (neg ok for geometric)
	plot.parms(parms.univ,c(0,5))

	parms <- parms.har; F <- parms$F; Finv <- parms$Finv; f <- parms$f;
	F(x)                  # wealth at x
	Finv(F(x)) - x        # should be 0

	f(x)                  # amount spent
	F(x) + f(x)           # wealth after spending
	
	Finv(F(x) + f(x))     # next position if do not reject
	Finv(F(x) + 0.5)      # new position if reject

	rm(Finv, F, f, parms)
	}
	
#####################################################################

simulate <- function(nSteps=5000) {
	# --- conditions
	p.reject <<- 0.01; 

	# --- iteration functions
	next.geo  <- build.recursion(parms.geo , p.reject, omega())
	next.har  <- build.recursion(parms.har , p.reject, omega())
	next.univ <- build.recursion(parms.univ, p.reject, omega())

	# --- initial values
	x0.geo  <- parms.geo$Finv(0.5); 
	x0.har  <- parms.har$Finv(0.5); 
	x0.univ <- parms.univ$Finv(0.5); 

	x.geo  <- rep(0,nSteps) ; x.geo [1]<-x0.geo
	x.har  <- rep(0,nSteps) ; x.har [1]<-x0.har
	x.univ <- rep(0,nSteps) ; x.univ[1]<-x0.univ

	for(i in 2:nSteps)	{ 
		x.geo [i] <- next.geo (x.geo [i-1]);
		x.har [i] <- next.har (x.har [i-1]);
		x.univ[i] <- next.univ(x.univ[i-1]);
		}
		z <- cbind(x.geo, x.univ, x.har)
		colnames(z) <- c("geo", "univ", "har")
		z 
}

# simres <- simulate()

# --- where does it stay
draw.hist <- function(x, parms) {
	hist(log(-parms$f(x[x<1000])), # avoid burn in
		breaks=100, main=parms$name, xlab="Alpha Levels",
		sub=paste("p=",p.reject,"  omega=",omega(),
		     "  mean(-f(x)) = ", round(-mean(parms$f(x)),5), "  sd=",round(sd(parms$f(x)),5))
		)
	}

see.results <- function() {
	
	# sequences
	i<- 1:500; parms <- parms.geo; idx <- "univ"
	
		plot(i, parms$F(simres[i,idx]), main=paste("Wealth(",idx,")"), type="p")
	 	plot(i,-parms$f(simres[i,idx]), main=paste("Alpha(",idx,")"), log="y")
	 

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




