
n <- 10000

r <- 0.5

x1 <- rnorm(n)
x2 <- r * x1 + sqrt(1-r^2)*rnorm(n)
y  <- (x1+x2+rnorm(n))/sqrt(3)

b<-(x1<1) & (x2<1); sum(b)

pairs(cbind(x1,x2,y),col=1+b)


out <- y[b]

qqnorm(out); abline(a=mean(out),b=sd(out))

mean(out)
sd(out)
