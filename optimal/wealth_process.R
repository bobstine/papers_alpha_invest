### Properties of the alpha investing processes

# --- recursion...get next position from prior
build.recursion <- function(parms, p, omega) { # p is rejection probability
	F <- parms$F; Finv <- parms$Finv; f <- parms$f; 
	function (x) { Finv(F(x) + ifelse(rbinom(1,1,p), omega, f(x))) }
	}
	

	
###  Harmonic  ###

K.har <- log(10000);

parms.har <- list(
	name = "Harmonic",
	F=      function(x) { K.har - log(x)}, 
	Finv =  function(F) { exp(K.har-F) },
	f =     function(x) { -1/x }
)


###  Universal ### (note you can solve for x < 1 to handle high wealth)

K.univ <- 1/2.10974

parms.univ <- list(
	name = "Universal",
	F = 	function(x) { K.univ/log(x+1) },
	Finv = function(F) { exp(K.univ/F)-1 },
	f = 	function(x) { -K.univ/((x+1)*log(x+1)^2) }
)


###  Geometric ### (note you can solve for x < 1 to handle high wealth)

psi <- 0.1
parms.geo <- list(
	name = "Geometric",
	F = 	function(x) { -psi*(1-psi)^(x-1)/log(1-psi) },
	Finv = function(F) { z <- log(1-psi); log(F*(psi-1)*z/psi)/z },
	f = 	function(x) { -psi*(1-psi)^(x-1) }
)

#####################################################################
# All have the problem of |f| > F for x < 1
# which is near where bid is equal to omega
# So, need to bound bid by omega

parms <- parms.univ

plot(parms$F,0,4); 
plot(function(x){-f(x)},0,20, add=TRUE, col="red")
f(1)

rm(parms)

#####################################################################
# --- check inverse functions
with( parms.har,  	Finv(F(10)))
with( parms.har,  	F(Finv(0.3)))

# --- plot to see wealth function F and bids (neg ok for geometric)
with( parms.univ,  plot(F,0.1,10))
with( parms.univ,  plot(f,0.1,100))


###### which(is.nan(x.univ))[1]
# 1.231302

Finv <- parms.univ$Finv; F <- parms.univ$F; f <- parms.univ$f … 

x <- 1.580502         # start
x <- 1.863090e+02
x <- 1.231302
x <- 0.544381  … At this point, start to spend more than you have… and it dies

F(x)                  # wealth at x
Finv(F(x)) - x        # should be 0

f(x)                  # amount spent
F(x) + f(x)           # wealth after spending
Finv(F(x) + f(x))     # next position if do not reject


Finv(F(x) + 0.5)      # new position if reject




rm(Finv, F, f)
#####
	
#####################################################################


# --- conditions
nSteps <- 10000; p <- 0.01; omega <- 0.5

# --- iteration functions
next.geo  <- build.recursion(parms.geo ,p,omega)
next.har  <- build.recursion(parms.har ,p,omega)
next.univ <- build.recursion(parms.univ,p,omega)

# --- initial values
x0.geo  <- parms.geo$Finv(0.5); 
x0.har  <- parms.har$Finv(0.5); 
x0.univ <- parms.univ$Finv(0.5); 

x.geo  <- rep(0,nSteps) ; x.geo [1]<-x0.geo
x.har  <- rep(0,nSteps) ; x.har [1]<-x0.har
x.univ <- rep(0,nSteps) ; x.univ[1]<-x0.univ

# --- simulate
for(i in 2:nSteps)	{ 
	x.geo [i] <- next.geo (x.geo [i-1]);
	x.har [i] <- next.har (x.har [i-1]);
	x.univ[i] <- next.univ(x.univ[i-1]);
	}


# --- where does it stay
draw.hist <- function(x, parms) {
	hist(x[x<1000], breaks=100, main=parms$name, xlab="Wealth Position",
		sub=paste("p=",p,"  omega=",omega,
		     "  f(med x)=", round(parms$f(median(x)),5), "  mean(f(x)) = ", round(mean(parms$f(x)),5))
		)
	}
	
par(mfrow=c(3,1))
	draw.hist(x.geo , parms.geo );
	draw.hist(x.har , parms.har );
	draw.hist(x.univ, parms.univ);
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





