### Properties of the alpha investing processes

###  Harmonic

K.har <- log(10000);

F.har    <- function(x) { K.har - log(x)}
Finv.har <- function(F) { exp(K.har-F) }
f.har    <- function(x) { -1/x }

###  Geometric  (note you can solve for x < 1 to handle high wealth)

psi <- 0.1
F.geo    <- function(x) { -psi*(1-psi)^(x-1)/log(1-psi) }
f.geo    <- function(x) { -psi*(1-psi)^(x-1) }
Finv.geo <- function(F) { z <- log(1-psi); log(F*(psi-1)*z/psi)/z }


# --- setup
name <- "Harmonic" ; F <- F.har; f <- f.har; Finv <- Finv.har
name <- "Geometric"; F <- F.geo; f <- f.geo; Finv <- Finv.geo

# --- check inverse functions
Finv(F(Finv(F(10))))

# --- plot to see wealth function F (neg ok for geometric)
plot(F,-2,100)

# --- get next position from prior
get.next <- function (x, p, omega) { # p is rejection probability
	return (Finv(
		F(x) + ifelse(rbinom(1,1,p), omega, f(x)) )
	); }
	
get.next(Finv(0.5),.01,.5)

hist.log <- function(h, main) {
	plot(h$mids, h$counts, log="x", pch=20, col="blue", type="h") }
	
#####################################################################

# --- set initial X location
x0 <- Finv(0.5); x0

# --- conditions
nSteps <- 20000; p <- 0.01; omega <- 0.5

# --- iterate
x <- rep(0,nSteps) ; x[1]<-x0
for(i in 2:nSteps)	x[i] <- get.next(x[i-1], p, omega)


# --- where does it stay
hist(x[x<1000], breaks=100, main=name); c(p, omega, median(x), f(median(x)), mean(f(x)))

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





