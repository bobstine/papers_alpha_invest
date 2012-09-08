###
###  p-value of normal test as function of mean
###


# --- distribution of the p-value under different alternative means

# returns p-values as function of percentile of distribution under Ha
# so the CDF of the p-values requires one to 'transpose' the figure
# Use your imagination!


# one-sided, reject if large mean
mu.a <- 0
plot( function(pct) (1-pnorm( qnorm(1-pct,mean=mu.a) )), .001, .999, log="y",
         ylim=c(.00000001,.5), xlab="Percentile", ylab="p-value"  )
for(m in 1:6) plot( function(pct) (1-pnorm( qnorm(1-pct,mean=m) )), .01, .99, add=TRUE )
   
# --- so if we set alpha = .0001, we get about 25% of the p-values (mu.a=3)
abline(h=0.0001,col="blue", lty=3)

#     whereas if they have .001 they get about 50%
abline(h=0.001,col="red", lty=3)
