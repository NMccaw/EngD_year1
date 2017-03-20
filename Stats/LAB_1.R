## MATH6157 Applied Statistical Modelling
## Introduction

## Slide 8, "A common example"
par(mai = c(1.36, 1.5, 1.093, 0.56))
curve(dnorm, -4, 4, ylab = expression(f[X](x)), cex.lab = 2, cex.axis = 2, lwd= 2) # plot a normal pdf
title("Probability density function", cex.main = 2) # add a title ...
legend("topleft", legend = c(expression(mu == 0), expression(sigma^2 == 1)), cex = 2) # ... and a legend
par(mai = c(1.36, 1.5, 1.093, 0.56))
curve(pnorm, -4, 4, ylab = expression(F[X](x)), cex.lab = 2, cex.axis = 2, lwd= 2) # plot a normal cdf
title("Cummulative distribution function", cex.main = 2)
legend("topleft", legend = c(expression(mu == 0), expression(sigma^2 == 1)), cex = 2)


## Slide 13, "Using R " 2
## ISwR by Dalgaard, p.110
library(ISwR) # load the package 'ISwR'
attach(thuesen) # attach a data.frame, for easy access to variables
## The thuesen data frame contains data on ventricular shortening velocity and blood glucose for 24 type 1 diabetic patients
velo.lm <- lm(short.velocity~blood.glucose, na.action = na.exclude) # fit a linear regression model
par(mai = c(1.36, 1.5, 1.093, 0.56))
plot(blood.glucose, short.velocity, cex.lab = 2, cex.axis = 2, cex = 2) # plot the data
lines(blood.glucose, fitted(velo.lm), lwd = 2) # add a regression line
segments(blood.glucose, fitted(velo.lm),
         blood.glucose, short.velocity, lty = 2, lwd = 2) # add vertical line segments showing the residuals
detach(thuesen) # detach the data.frame

## Slide 20, "Descriptive statistics"
y <- rnorm(50) # generate some random data from a standard normal distribution
mean(y) 
sd(y)
var(y)
median(y)

## Slide 21, "Qunatile and five-number summary"
summary(y) # five-number summary (min, Q1, median, Q2, max) along with the mean
pvec <- seq(0, 1, 0.1) # sequence from 0 to 1 in steps of 0.1
pvec
round(quantile(y, pvec), 3) # quantiles using the percentage points from pvec (deciles), rounded to 3 significant figures

## Slide 22, "Summaries by group"
library(ISwR) # load the package 'ISwR'
attach(energy) # attach the energy data.frame
## the energy data frame contains data on energy expenditure during exercise by two groups of lean and obese women.
tapply(expend, stature, mean) # apply the function mean to the variable 'expend', grouping by the levels of 'statue'
tapply(expend, stature, length) # ditto for length, to calulate the number in each group
tapply(expend, stature, sd) # and again, to calculate the standard deviation in each group
## with better formatting
xbar <- tapply(expend, stature, mean)
s <- tapply(expend, stature, sd)
n <- tapply(expend, stature, length)
cbind(mean = xbar, std.dev = s, n = n) # cbind binds vectors together into an array
detach(energy) # detach the data.frame

## Slide 23, "Graphical visualisation of data"
par(mai = c(1.36, 1.5, 1.093, 0.56))
hist(y, breaks = 6, cex.main = 2, cex.lab = 2, cex.axis = 2, col = "blue") # histogram, with 6 bins

## Slide 24, "Empirical cumulative distribution"
n <- length(y) # how many elements in y?
par(mai = c(1.36, 1.5, 1.093, 0.56))
plot(sort(y), (1:n)/n, type = "s", ylim = c(0, 1), lwd = 2, cex.lab = 2, cex.axis = 2,
     xlab = "Sorted data", ylab = "ecdf") # plotting the ecdf
## sort sorts into ascending order by default; 1:n produces a sequence from 1 to n in steps of size 1; 
## c(0, 1) is a vector with entries 0 and 1

## Slide 25, "Graphs of grouped data"
library(ISwR)
attach(energy)
expend.lean <- expend[stature == "lean"] # subset the data in expend to only include those observations in the 'lean' group
expend.obese <- expend[stature == "obese"] # as above, but subsetting to only include observations in the 'obese' group
par(mfrow=c(2, 1)) # produce two plots, arranged as two rows and one column
par(mai = c(1.36, 1.5, 1.093, 0.56))
hist(expend.lean, breaks = 10, xlim = c(5, 13), ylim = c(0, 4), col = "white", cex.main = 2, cex.lab = 2, cex.axis = 2)
hist(expend.obese, breaks = 10, xlim = c(5, 13), ylim = c(0, 4), col = "grey", cex.main = 2, cex.lab = 2, cex.axis = 2)
par(mfrow=c(1, 1)) # reset the plotting device
detach(energy)

## Slide 27, "Graphs of grouped data"
attach(energy)
par(mai = c(1.36, 1.5, 1.093, 0.56))
boxplot(expend ~ stature, ylab = "expend", cex.lab = 2, cex.axis = 2, lwd = 2) # boxplots of 'expend' for each level of 'stature'
detach(energy)

## Slide 28, "Higher dimensional data"
attach(thuesen)
par(mai = c(1.36, 1.5, 1.093, 0.56))
plot(blood.glucose, short.velocity, pch = 16, cex.lab = 2, cex.axis = 2, cex = 2) # plot the data, pch controls the plotting character
detach(thuesen) # detach the data.frame

## Slide 29, "Higher dimensional data"
par(mai = c(1.36, 1.5, 1.093, 0.56))
pairs(iris[1:4], main = "Anderson's Iris Data -- 3 species",
      pch = 16, col = c("red", "green3", "blue")[iris$Species], cex.lab = 2, cex.axis = 2, cex.main = 2)
## 'pairs' produces a scatter-plot matrix, with 2D scatter plots between each pair of variables.
## We have also colour coded the points by 'species'

## Slide 35/36, "Empirical CDFs"
y <- rnorm(50) # generate some random data
y.ecdf <- ecdf(y) # estimate the cdf from the data
par(mai = c(1.36, 1.5, 1.093, 0.56))
plot(y.ecdf, xlab = "y", ylab = "F(y)", main = "", cex.lab = 2, cex.axis = 2, lwd = 2) # plot the ecdf
a <- 0 
y.ecdf(a) # P(Y < a)
a.ecdf <- y.ecdf(0)
lines(c(a, a), c(-10, a.ecdf), lty = 2, lwd = 2) # display this probability graphically
lines(c(-10, a), c(a.ecdf, a.ecdf), lty = 2, lwd = 2)

y2 <- rbinom(50, 10, .6)
y2.ecdf <- ecdf(y2)
par(mai = c(1.36, 1.5, 1.093, 0.56))
plot(y2.ecdf, xlab = "y", ylab = "F(y)", main = "", cex.lab = 2, cex.axis = 2, lwd = 2) # plot the ecdf
a2 <- 7 
y2.ecdf(a2) # P(Y < a)
a2.ecdf <- y2.ecdf(a2)
lines(c(a2, a2), c(-10, a2.ecdf), lty = 2, lwd = 2) 
lines(c(-10, a2), c(a2.ecdf, a2.ecdf), lty = 2, lwd = 2)

## Slide 38/39 "Two (more) examples"
y <- 0:50 # sequence from 0 to 50
par(mai = c(1.36, 1.5, 1.093, 0.56))
plot(y, dbinom(y, size = 50, prob = .33), type = "h", ylab = "p(y)", cex.axis = 2, cex.lab = 2, lwd = 2) # plot the binomial probabilities with m=50 and rho=.33

par(mai = c(1.36, 1.5, 1.093, 0.56))
curve(dexp, 0, 10, ylab="f(y)", xlab = "y", lwd = 2, cex.axis = 2, cex.lab = 2)

## Slide 43 "Distributions in R"
mu <- 5
dpois(7, mu) # probability mass function (probability Y = 7)
ppois(7, mu) # cumulative distribution function (probability Y <= 7)
qpois(.5, mu) # returns the 0.5-quantile (median) - inverts the cdf
rpois(10, mu) # generate 10 (pseudo) random numbers from Poisson(10)

## Slide 52 "Examples"
## functions, which both calculate the **minus** log-likelihood for the binomial distribution
## Note the use of vectorisation in both functions
## - lchoose, log and dbinom are all applied to vectors (e.g. y) and so no explicit looping is required
mloglik1 <- function(prob, y, m) { # "from scratch"
  l <- lchoose(m, y) + y * log(prob) + (m - y) * log(1 - prob)
  -1 * sum(l)
}

mloglik2 <- function(prob, y, m){ # uses the dbinom function to calculate the binomial probabilities
  -1 * sum(log(dbinom(y, m, prob)))
}

n <- 20 # size of data set
m <- 10 # number of binomial trials
tp <- 0.7 # true probability
y1 <- rbinom(n, m, tp) # generate some binomial data
mle1 <- nlm(m?nlmloglik2, p = .5, y = y1, m = m)
## you will probably get warnings; the binomial log-likelihood function isn't always well behaved.
sum(y1)/(n * m) # true MLE

## Slide 53 "Examples"
## comparing the log-likelihood for different sample sizes
y2 <- rbinom(5 * n, m, tp) # generate some (more) binomial data (bigger sample size)
mle2 <- nlm(mloglik2, p = .5, y = y2, m = m) 
p <- seq(.05, .95, by = .01) # sequence of probabilities
l1 <- sapply(p, mloglik2, y = y1, m = m) # sapply allows us to apply the function "loglik2" to each element of p
l2 <- sapply(p, mloglik2, y = y2, m = m)
par(mai = c(1.36, 1.5, 1.093, 0.56))
matplot(p, cbind(-1 * l1 + mle1$minimum, -1 * l2 + mle2$minimum), type = "l", ylab = "Log-likelihood", lwd = 2, cex.axis = 2, cex.lab = 2, 
        col = c("black", "red")) # plot the results, standardised by subtracting the maximium log-likelihood
abline(v = sum(y1)/(n * m), lwd = 2, lty = 3) # add in vertical lines at the actual MLEs
abline(v = sum(y2)/(5 * n * m), lwd = 2, lty = 4, col = "red")

## Slide 55 "Examples", normal distribution
mloglik3 <- function(mu_sigma2, y) { ## calculate the log-likelihood; first argument is vector c(mu, sigma2)
  mu <- mu_sigma2[1]
  sigma <- sqrt(mu_sigma2[2])
  -1 * sum(dnorm(y, mu, sigma, log = T))
}
n <- 100
y3 <- rnorm(n, mean = 1, sd = 2)
mle3 <- nlm(mloglik3, p = c(0, 4), y = y3)
mle3
mean(y3)
(n-1) * var(y3)/n

## Slide 57 "Asymptotic distribution of MLEs"
## back to the binimial example
mle1 <- nlm(mloglik2, p = .5, y = y1, m = m, hessian = T) # rerun the minimisation but now return the hessian
phat <- mle1$estimate # MLE
phat
ofi <- mle1$hessian # observed Fisher information (we don't need to multiple by -1, as we are working with the minus log-lik.)
ofi
conf.level <- 0.95
crit <- qnorm((1+conf.level)/2) # pick the correct quantile of the normal distribution
phat + c(-1, 1) * crit * sqrt(1/ofi) # confidence interval

## Slide 59, hypothesis testing
## using our previous normal data
mloglik.null <- function(sigma2, y) { ## calculate the log-likelihood under H0: mu=0; first argument is vector sigma2
  -1 * sum(dnorm(y, 0, sqrt(sigma2), log = T))
}
mle.null <- nlm(mloglik.null, p = 1, y = y3) # find mle for sigma2 under H0
mle.alt <- mle3 # mles for mu and sigma2 under H1
LR <- - 2 * ( -1 * mle.null$minimum + mle.alt$minimum) # don't forget that our loglikelihood functions return *-loglik*!
LR
1 - pchisq(LR, 1) # how likely is this ratio under the reference distribution?



