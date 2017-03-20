library(faraway)
data("pima")
pima
# use summary to view quick data to see if anything is out of place
summary(pima)

sort(pima$diastolic)

#set all data with values 0 to NA
pima$diastolic[pima$diastolic == 0] <- NA
pima$glucose[pima$glucose == 0] <- NA
pima$triceps[pima$triceps == 0] <- NA
pima$insulin[pima$insulin == 0] <- NA
pima$bmi[pima$bmi == 0] <- NA
# use factor to treat the variable as categorical not quantatative
pima$test <- factor(pima$test)

summary(pima$test)
# use levels to use descriptive labels
levels(pima$test) <- c("negative","positive")

hist(pima$diastolic)
# histograms can be misleading. therefore use kernal density functions which are essentially 
# smoothed histomgrams
plot(density(pima$diastolic,na.rm=TRUE))
# can also plot the sorted data
plot(sort(pima$diastolic),pch=".")
#here are two bivariable plots

#1st a standard scatterplot
plot(diabetes ~ diastolic,pima)
# 2nd boxplot suitable for showing a quantitative and qualitative variables
plot(diabetes ~ test,pima)
# can also have a scatter plot matrix using
pairs(pima)

# Regression

#look at relationship between final scores and mid term scores
data("stat500")
stat500 <- data.frame(scale(stat500))
plot(final ~ midterm,stat500)
abline(0,1)
# add least square fit line using lm
g <- lm(final ~ midterm,stat500)
abline(g$coef,lty=5)
#also can compute correlations
cor(stat500)

## Residual diagnostics

attach(faithful)
yellowstone<-lm(eruptions~waiting)
plot(waiting,eruptions)
abline(yellowstone$coef)


# Exercise 1: tomatoes

Fertelizer = c(5,6,8,11,12,13,14,15,16,17,18)
Yield = c(18,20,21,25,24,27,30,31,31,33,29)
tomatoe = lm(Yield~Fertelizer)
plot(Fertelizer,Yield)
abline(tomatoe$coef)
cor(Yield,Fertelizer)
summary(tomatoe)
plot(tomatoe)

# Exercise 2: Lactic Acid

x = c(1,1,1,1,3,3,3,3,5,5,5,5,10,10,10,10,15,15,15,15)
y = c(1.1,0.7,1.8,0.4,3.0,2.4,4.9,4.4,4.5,7.3,8.2,6.2,12.0,13.1,12.6,13.2,18.7,19.7,17.4,17.1)
plot(x,y)
model = lm(y~x)
lactic.lm<-lm(y~ x-1)
abline(lactic.lm$coefficients)
plot(lactic.lm)
plot(resid(lactic.lm))

sumary(lactic.lm)