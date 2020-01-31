### As discussed, it seems like R may require these to be in a slightly different form
### While we are working on that, I wanted to confirm the functions were working.

R = 1.987  ### Gas Constant

##############################
### Single State Functions ###
##############################

SingleStateMeltingFunction <- function(T, Sn, Su, dH, Tm)
{
  (Sn + Su*exp(-(dH/R)*(1/Tm - 1/T)))/(1+exp(-(dH/R)*(1/Tm - 1/T)))
}

SingleStateFractionSn <- function(T,dH,Tm)
{
  1/(1 + exp(-((dH/R)*((1/Tm) - (1/T)))))
}

SingleleStateFractionSu <- function(T,dH,Tm)
{
  (exp(-((dH/R)*((1/Tm) - (1/T)))))/(1 + exp(-((dH/R)*((1/Tm) - (1/T)))))
}

### nonlinear fit to single state

raw.data <- read.table('2GKU_Melt_2_CD.txt')
temperatures <- seq(4, 94, 2) + 273.15
wavelengths <- seq(220, 340, 1)

plot(wavelengths, raw.data[,1])
plot(temperatures, raw.data[15,])

svd.temp <- svd(raw.data)

plot(wavelengths, -svd.temp$u[,1])
plot(wavelengths, -svd.temp$u[,2])
plot(wavelengths, -svd.temp$u[,3])

plot(temperatures, svd.temp$v[,1]*svd.temp$d[1])
plot(temperatures, svd.temp$v[,2])
plot(temperatures, svd.temp$v[,3])

y <- svd.temp$v[,1]*svd.temp$d[1]
x <- temperatures
model.onestate <- y ~ (Sn + Su*exp(-(dH/R)*(1/Tm - 1/x)))/(1+exp(-(dH/R)*(1/Tm - 1/x)))
m<-nls(model.onestate,
    start = list(Sn = -100, Su = -1, dH = -30000, Tm=325))
summary(m)
predict(m, x)
plot(temperatures, y)
lines(temperatures, predict(m,x))

y <- svd.temp$v[,2]*svd.temp$d[2]
x <- temperatures
m<-nls(y ~ (Sn + Su*exp(-(dH/R)*(1/Tm - 1/x)))/(1+exp(-(dH/R)*(1/Tm - 1/x))),
       start = list(Sn = -100, Su = -1, dH = -30000, Tm=325))
summary(m)
predict(m, x)
plot(temperatures, y)
lines(temperatures, predict(m,x))

y <- svd.temp$v[,3]*svd.temp$d[3]
x <- temperatures
m<-nls(y ~ (Sn + Su*exp(-(dH/R)*(1/Tm - 1/x)))/(1+exp(-(dH/R)*(1/Tm - 1/x))),
       start = list(Sn = 5, Su = -15, dH = -20000, Tm=325))
summary(m)
predict(m, x)
plot(temperatures, y)
lines(temperatures, predict(m,x))

### BIG QUESTION : How can we fit all 3 of these at one time.  Let dH, Tm be fixed, but different S* for each.

##############################
### Two State Functions ###
##############################

TwoStateMeltingFunction <- function(T,Sn,S1,Su,dH1,Tm1,dH2,Tm2) #Native to Intermediate 1 to Unfolded
{
  (Sn + S1*exp(-((dH1/R)*((1/Tm1) - (1/T)))) + 
     Su*exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T)))) +
                                                                                exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))
  
}

model.twostate <- y ~ TwoStateMeltingFunction(T, Sn, S1, Su, dH1, Tm1, dH2, Tm2)

### nonlinear fit to two states

raw.data <- read.table('2GKU_Melt_2_CD.txt')
temperatures <- seq(4, 94, 2) + 273.15
wavelengths <- seq(220, 340, 1)

plot(wavelengths, raw.data[,1])
plot(temperatures, raw.data[15,])

svd.temp <- svd(raw.data)

plot(wavelengths, -svd.temp$u[,1])
plot(wavelengths, -svd.temp$u[,2])
plot(wavelengths, -svd.temp$u[,3])

plot(temperatures, svd.temp$v[,1]*svd.temp$d[1])
plot(temperatures, svd.temp$v[,2])
plot(temperatures, svd.temp$v[,3])

y <- svd.temp$v[,2]*svd.temp$d[2]
T <- temperatures
plot(T,y)

m<-nls(model.onestate,
       start = list(Sn = -5, S1 = 25, Su = 10, dH1 = -20000, Tm1=320, dH2 = -20000, Tm2 = 335))
summary(m)

plot(temperatures, y)
lines(temperatures, predict(m,x))



#### nlm 

f <- function(x) sum((x-1:length(x))^2)
nlm(f, c(10,10))
nlm(f, c(10,10), print.level = 2)
utils::str(nlm(f, c(5), hessian = TRUE))


f <- function(T)
