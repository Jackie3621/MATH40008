#read in the data
df <- read.csv("GOOG.csv", header = TRUE, sep = ',')
#extract the required column
S <- df$Close
#define dt as 1/no. of trading days per year
dt <- 1/252
#calculate the log return
R <- diff(log(S))

#estimation from empirical data
muhat <- (mean(R)/dt+(var(R))/(2*dt))
sigmahat <- sqrt(var(R)/dt)

#for BS model
mu <- (muhat-0.5*(sigmahat)^2)*dt
sigma <- sigmahat*sqrt(dt)

#threshold for a jump
upper2 <- mu + 2*sigma
lower2 <- mu - 2*sigma

par(mar = c(5, 5, 1, 1))

#Figure 1
plot(S, type = 'l', ylab = 'price',bty='L',xlab='day', xlim= c(0,504), ylim=c(0,3100),
     xaxs="i", yaxs="i",cex.lab = 1.5,cex.axis = 1.5)

#Figure 2
plot(R, type = 'h', xlab = 'day', ylab = 'empircial log-returns',
     xaxs= 'i',ylim=c(-0.12,0.1),cex.lab = 1.5,cex.axis = 1.5)
abline(upper2, 0, col = 'blue')
abline(lower2, 0, col = 'blue')

#set values within thresholds to be zero
r <- R
r[upper2 > r & lower2 < r] <- 0
#index of date which jump happens
date <- which(r != 0)

#MLE for HPP
#estimation analytically of theta = 1/mean(X)
T <- diff(c(0, date, 504))
#log likelihood function
log_like <- function(theta, X){
  n <- length(X)
  loglik <- n*log(theta) - n*theta*mean(X)
  return(-loglik)
}
MLE <- optim(par = c(0.1), log_like, X = T, method = "Brent", lower = 0, upper = 100)

#Figure 3
pdf(file = "R.pdf", width = 6, height = 4)
par(mar = c(5, 5, 1, 1))
plot(c(0,1:length(date),length(date))~c(0,date,504), type = 's', col='blue',xlab = 'day',
     ylab = 'Poisson process N(t)',cex.lab = 1.5,cex.axis = 1.5)
legend("bottomright",inset=.05, legend=c("empircial data", "simulations"),
       col=c("blue", "black"), lty=c (1,1), cex=1.1)
set.seed(34)
B <- 504
lambda <- MLE$par
m <- round(3*B*lambda)
u1 <- runif(m,0,1)
u2 <- runif(m,0,1)
u3 <- runif(m,0,1)
t1 <- -1/lambda*log(u1)
t2 <- -1/lambda*log(u2)
t3 <- -1/lambda*log(u3)
s1 <- cumsum(t1)
s1 <- s1[s1<B]
s2 <- cumsum(t2)
s2 <- s2[s2<B]
s3 <- cumsum(t3)
s3 <- s3[s3<B]
n1 <- length(s1)
n2 <- length(s2)
n3 <- length(s3)
Nt1 <- 1:n1
Nt2 <- 1:n2
Nt3 <- 1:n3
lines(c(0,Nt1,n1)~c(0,s1,B), type = 's')
lines(c(0,Nt2,n2)~c(0,s2,B), type = 's')
lines(c(0,Nt3,n3)~c(0,s3,B), type = 's')
u4 <- runif(m,0,1)
u5 <- runif(m,0,1)
u6 <- runif(m,0,1)
t4 <- -1/lambda*log(u4)
t5 <- -1/lambda*log(u5)
t6 <- -1/lambda*log(u6)
s4 <- cumsum(t4)
s4 <- s4[s4<B]
s5 <- cumsum(t5)
s5 <- s5[s5<B]
s6 <- cumsum(t6)
s6 <- s6[s6<B]
n4 <- length(s4)
n5 <- length(s5)
n6 <- length(s6)
Nt4 <- 1:n4
Nt5 <- 1:n5
Nt6 <- 1:n6
lines(c(0,Nt4,n4)~c(0,s4,B), type = 's')
lines(c(0,Nt5,n5)~c(0,s5,B), type = 's')
lines(c(0,Nt6,n6)~c(0,s6,B), type = 's')
dev.off()

#Figure 4
qqplot(x=qexp(ppoints(length(T))), y=T, 
     xlab="Theoretical Quantiles", ylab= "Empircial Quantiles",cex.lab = 1.5,cex.axis = 1.5)
qqline(T, distribution=qexp, col='blue')

#KS test to reassure QQ plot
ks.test(T, "pexp", rate=MLE$par)

#MLE for Hawkes
test <- function(x, data = t){
  if(x[3] <= x[2]){
    return(Inf)
  }
  r <- likelihoodHawkes(x[1], x[2], x[3], data)
  return(r)
}
MLE3 <- optim(c(MLE$par,1,2), test, data = date)

#Figure 5
dS <- abs(replace(diff(S), c(1:503)[-date],0))
size <- replace(dS,date,MLE3$par[2])
inten <- function(t){
  t_k <- c(1:t)
  lambda_t <- MLE3$par[1] + sum(size[1:t]*exp(-MLE3$par[3]*(t-t_k)))
  return(lambda_t) 
}

intensity <- c()
for(i in 1:504){
  intensity[i] <- inten(i)
}

plot(intensity,type='l',ylab = 'estimated intensity',bty='L',xlab='day', xlim= c(0,504),
     xaxs="i", col = 'blue',cex.lab = 1.5,cex.axis = 1.5)
abline(MLE$par,0, col = 'red')
abline(MLE3$par[1],0,lty = 3)
legend("topright",inset=.02, legend=c("HPP", "Hawkes process", "reversion level of Hawkes process"),
       col=c("red", "blue", "black"), lty=c (1,1,2), cex=1.1)

#AIC = 2k - 2ln(L)
AIC1 <- 2*1 + 2*MLE$value
AIC3 <- 2*3 + 2*MLE3$value

#MLE3
test(MLE3$par,date)
#2018-2019
test(c(0.03869434,0.04351336,0.10704307),date)
#2015-2019
test(c(0.02797943,0.04713436,0.1054357),date)
#2010-2019
test(c(0.02865434,0.03691492,0.10352264),date)
#2005-2019
test(c(0.01623701,0.04530607,0.06725749),date)