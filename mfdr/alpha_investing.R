# $Id: alpha_investing.R,v 1.17 2007/07/23 18:13:49 bob Exp $
# R code for alpha-investing

#
#  NEED: add an upper bound for the amout of alpha-wealth can spend (say 1)
#

# An alpha-investing procedure combines an investing rule with a regulator
# function.  In order to have local access to the test results, these functions
# have to be defined within the scope of the make.test.procedure function.
#
# Arguments to make.test.procedure are W(0), \omega, and the investing rule
# The arguments to the rule are the current wealth and the state of the procedure.
# The state of the test procedure is the triple (current.index, number since rejected, 
# and the number rejected).  Eta is set to 0.95 globally by default.
#
# If m, the number of tests is specified as 0, then the procedure assumes it will
# be applied to a possibly infinite stream of hypotheses, none of which are revisited.
#
# Regulators take as input alpha level, p-value, current wealth, and the known lower
# bound for p-value, which is zero in most cases unless the H0 has been test previously
# as in the 'nibbling' approach.
#

.eta <- 0.95

threshold.from.alpha <- function(alpha, lower.bound) {
	lower.bound + alpha * (1-lower.bound)
}
			
log.regulator <- function(al, p, omega, lower.bound=0) { 
  	threshold <- threshold.from.alpha(al, lower.bound);
  	if(p <= threshold) {
  		level <- (p-lower.bound)/(1-lower.bound); # need 'conditional' p-value
  		omega + log(1-level) }
  	else log(1-al) 
}

recip.regulator <- function(al, p, omega, lower.bound=0) {
  	threshold <- threshold.from.alpha(al, lower.bound);
  	if(p <= threshold)  omega 
  	else -al/(1-al) 
}

default.regulator <- recip.regulator;

make.test.procedure <- function(m, alpha.function, 
                       regulator = default.regulator, wealth=0.05*.eta, payout=0.05) {
	r           <- vector(mode="logical", length=m);   # rejections
	lower.bound <- vector(mode="numeric", length=m);   # condition on p-value larger than this
	state <- list (
		current.index = 0,      	# jth test
		number.since.rejected = 0,	# tests since last rejected
		number.rejected = 0);		# cumulate number rejected
	update.state <- function(threshold, p, d) {	
		# cat("Test at level alpha=",a," vs p=",p," with delta = ", d, "\n");
		wealth <<- wealth + d;
		if (p<=threshold) { # reject
			if(m > 0) r[state$current.index] <<- TRUE    # mark location
			else      r <<- c(r,state$current.index);    # store index
			state$number.rejected <<- state$number.rejected + 1;
			state$number.since.rejected <<- 0; }
		else { # not rejected
			if (m > 0) lower.bound[state$current.index] <<- threshold;
			state$number.since.rejected <<- state$number.since.rejected + 1; }
	};
	list(
		test = function(j, pj) {
			state$current.index <<- j;
			al <- min(.5, max(0, alpha.function(wealth, state)));
			bound <- ifelse(m > 0, lower.bound[j], 0);
			delta <- regulator(al, pj, payout, bound);
			while((wealth+delta) < -0.01) {   # make sure do not overspend
				cat("Patching... w+d = ",wealth,"+",delta,"=",wealth+delta,"\n");
				cat("al=",al," pj=",pj," payout=",payout," bound=",bound,"\n");
				al<- al/2; 
				delta <- regulator(al, pj, payout, bound) 
			}
			p.threshold <- threshold.from.alpha(al,bound);
			update.state(p.threshold, pj, delta);
			return(pj <= p.threshold)
		},
		wealth = function () { wealth },
		test.results = function () { r },
		number.rejected = function () { state$number.rejected },
		state = function () { state },
		print = function() {
			cat("State: j= ", state$current.index, 
			          " t= ", state$number.since.rejected,
			          " R= ", state$number.rejected, 
			          " W = ", wealth,
			          " a = ", alpha.function(wealth,state), "\n");
		}
	)
}


make.bonferroni.test.procedure <- function(m, 
                                 regulator=default.regulator, wealth=0.05*.eta, payout=0.05) {
	alpha.rule <- function(...) { wealth/m; }
	make.test.procedure(m, alpha.rule, regulator=regulator, wealth=wealth, payout=payout)
}

make.aggressive.test.procedure <- function(m, 
                                  regulator = default.regulator, wealth=0.05*.eta, payout=0.05) {
	if(m > 0) { # linear decay over remaining
		alpha.rule <- function(w,state) { # this version spends all by test m
			max(w/(2+state$number.since.rejected), w/(w+1+m-state$current.index)) }}
	else alpha.rule <- function(w,state) w/(2+state$number.since.rejected) 
	make.test.procedure(m, alpha.rule, regulator=regulator, wealth=wealth, payout=payout)
}

aggressive.sequential.test <- function(p, regulator=default.regulator, alpha=0.05) {
	m <- length(p)
	test <- make.aggressive.test.procedure(m, wealth=alpha*.eta, regulator=regulator);
	# one pass is all you need.
	for (j in 1:m){ test$test(j,p[j]); }  # add test$print() to see results inside
	test$test.results()
}

# show what goes wrong without the ordering; cost of aggressive
shuffled.aggressive.sequential.test <- function(p, regulator=default.regulator, alpha=0.05) {
	m <- length(p)
	# permute indices
	j <- sample(1:m,m,replace=FALSE)
	r <- aggressive.sequential.test(p[j], regulator=regulator, alpha=alpha)
	r[order(j)]
}



# Simes is the nibbling approach, with closer ties to BH method
simes.sequential.test <- function(p, regulator=default.regulator, alpha=0.05) {
	m <- length(p)
	alpha.level <- 0;  # need protected global variable in following function
	alpha.rule <- function(...) { alpha.level } 
	test <- make.test.procedure(m, alpha.rule, wealth=alpha*.eta, regulator=regulator);
	while( (w<-test$wealth()) > 0.001 && (nr<-test$number.rejected())<m ) {
		alpha.level <- w/(w+m-nr); #  a/(1-a) penalty
		# cat("Simes: level = ", alpha.level, " nr = ", nr, "  "); test$print();
		not.rejected <- (1:length(p))[!test$test.results()];
		for (j in not.rejected){ test$test(j,p[j]); }
	}
	test$test.results()
}


##############################################################################
#   Simple to program methods
##############################################################################
# 
# These tests return a vector indicating which tests were rejected
#

bh.test <- function(p, alpha=0.05) {
	m <- length(p);
	ord <- order(p);
	i <- 1;
	r <- vector("logical", m);
	al <- alpha/m;
	while(p[ord[i]] <= i*al)  {
		r[ord[i]] <- TRUE; i <- i + 1; 
	}
	r
}

# assumes false come first as produced by generate.pvalues
weighted.p.values <- function(p,m1) {
	m <- length(p);
	p[1:m1] <- p[1:m1]/(m/m1);
	p[(m1+1):m] <- 1;
	return(p)
	}
	
sim.wbh.test <- function(p, alpha=0.05) {
	# assume ordered p-values as generated by generate.pvalues
	return(bh.test(weighted.p.values(p,.number.false.nulls),alpha))
	}
 	
# false.null indicates which are false
wbh.test <- function(p, false.null, alpha=0.05) {
	m <- length(p);
	w <- rep(0.0000000001,m)
	# assign positive weight to those with false null, "zero" to others
	w[false.null] <- m/sum(false.null);
	# cat("sum w = ", mean(w),"\n");
	return(bh.test(p/w,alpha=alpha))
	}



bonferroni.test <- function(p, alpha=0.05) { p <= alpha/length(p) }

naive.test <- function(p, alpha=0.05) { p <= alpha }

##############################################################################
#
# Simulate tests: stores in *global* v and r
#
#    Test procedure has two arguments: vector of p-values, scalar alpha
# 
##############################################################################

generate.pvalues <- function(m, pi=0.5, sigma=sqrt(2*log(m))) {
	m1 <- floor(m * pi);
	.number.false.nulls <<- m1
	mu <- c(sort(abs(rnorm(m1, sd=sigma)), decreasing=TRUE),vector("integer",m-m1));
	1-pnorm(mu+rnorm(m))
	}
	
# p <- generate.pvalues(10, 5); p

simulate.tests <- function(tests, m, pi) {
	N <- nrow(v);
	n.tests <- length(tests);
	cat("Simulating ", n.tests," tests for ", N," reps. [m=",m,",pi=",pi,"]\n");
	true.null.position <- ((floor(m*pi))+1):m;
	for(rep in 1:N) {
		p <- generate.pvalues(m, pi=pi);
		for(i in 1:n.tests) {
			reject <- tests[[i]](p);
			r[rep,i] <<- sum(reject);
			v[rep,i] <<- sum(reject[true.null.position] != FALSE);
		}
	}
}


##############################################################################
#
# Simulate streaming tests
# 
##############################################################################

# simple function for fdr calculation
fdr <- function(v,r) { ifelse(r>0,v/r,0) }

# probability of transition from state 0 (null) to state 1 (mean mu)
run.experiment <- function(p01.list, odds, n.trials) {
	for(p01 in p01.list) {	
		config <- list(p01=p01, odds=odds, n.trials=n.trials);
		results <- simulate.configuration(config$p01, config$odds, config$n.trials);
		config$n       <- mean(results[,"n"]);  # number of tests done
		config$wealth  <- mean(results[,"W"]);  # wealth at end
		config$fdr     <- mean(results[,"FDR"]);
		config$power   <- mean((results[,"R"]-results[,"V"])/results[,"S1"]);
		config$fpower  <- mean(results[,"ffirst"]/results[,"nfirst"]);
		config$lpower  <- mean(results[,"flast" ]/results[,"nlast" ]);
		config$nfirst  <- mean(results[,"nfirst"]);
		config$nlast   <- mean(results[,"nlast"]);
		simres <<- c(simres, list(config));
	}}

simulate.configuration <- function(p01, odds, n.trials, mu1=3) {
	p10 <- odds * p01;
	p00 <-  1-p01; p11 <- 1-p10;   	 
	results <- matrix(nrow=n.trials, ncol=10)
	colnames(results) <- c("n","W", "FDR","V","R","S1", "nfirst","nlast","ffirst","flast")
	max.n.tests <- 4000;      	# room for intermediate results
	for(trial in 1:n.trials) {  # grab properties of test (eg, wealth) from the test
		test <- make.aggressive.test.procedure(0); #  make.test.procedure(0, alpha.rule); 
		total.V <- total.R <- total.State1 <- 0
		num.first <- num.last <- found.first <- found.last <- 0
		j <- 1; prior.state <- 0; state <- 1
		while(j <= max.n.tests && test$wealth()>1.0e-10) {
			next.state <- ifelse(state==1, rbinom(1,1,p11), 1-rbinom(1,1,p00))
			if (state > prior.state) num.first <- num.first+1
			if (state > next.state)  num.last  <- num.last +1
			total.State1 <- total.State1 + state;
			z <- rnorm(1) + ifelse(state==1,mu1,0);
			p <- 2*pnorm(-abs(z));
			reject <- test$test(j,p);
			if(reject) {  # cat("p = ", p, " in state ", state, " with z = ", z,"\n");
				if (state > prior.state) found.first <- found.first+1  # found first
				if (state > next.state)  found.last  <- found.last +1  # found last
				total.R <- total.R + 1;
				total.V <- total.V + 1 - state; # rejected but should not have
			}
			# debug s[j] <- state; if(reject) { r[j] <- 1; v[j] <- 1 - state; }
			prior.state <- state;	state <- next.state; 
			j <- j + 1;
		}
		results[trial,] <- c(test$state()$current.index, test$wealth(),
								fdr(total.V,total.R),total.V,total.R,total.State1,
								num.first, num.last, found.first, found.last)
	}
	return(results)
}


extract <- function(lst, component.name, match.list=NULL) {
	if(0==length(match.list))
		unlist(lapply(lst,function(x) x[[component.name]]))
	else {
		match.name <- names(match.list)[1];
		matches <- unlist(lapply(lst,function(x) x[[match.name]]==match.list[[match.name]]));
		if (1 == length(match.list)) extract(lst[matches], component.name)
		else extract(lst[matches], component.name,match.list[2:length(match.list)])
	}}


