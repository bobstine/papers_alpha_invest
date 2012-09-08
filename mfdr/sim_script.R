# R script for grid simulation pipe
source("alpha_investing.R");
simres <- vector("list", 0);
p01 <- c(.05,.025,.01,.005,.0025);
n.trials <- 1000;
odds <- 9;
run.experiment(p01, odds, n.trials);
save(simres, file=paste("simres_",odds,sep=""));
