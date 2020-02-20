### As discussed, it seems like R may require these to be in a slightly different form
### While we are working on that, I wanted to confirm the functions were working.


### Setup ####

### Data

raw.data <- read.table('2GKU_Melt_2_CD.txt')
temperatures <- seq(4, 94, 2) + 273.15 ## Temperature in Kelvin
wavelengths <- seq(220, 340, 1)

### Visualize, maybe a 3d plot?

plot(wavelengths, raw.data[,1])
plot(temperatures, raw.data[55,])

### Need to run 'estimated' parameters simulation as in Mathematica

### 
svd.temp <- svd(raw.data)

temp.2 <- svd.temp$u%*%diag(svd.temp$d)%*%t(svd.temp$v)

plot(temperatures, temp.2[55,])

singular.temp <- svd.temp$d
singular.temp

singular.temp[3:46]<-0

singular.temp

temp.3 <- svd.temp$u%*%diag(singular.temp)%*%t(svd.temp$v)

plot(temperatures, temp.2[55,], type='l')
plot(temperatures, temp.3[55,], type='l')

cbind(raw.data[,1],temp.2[,1], temp.3[,1])

plot(wavelengths, -svd.temp$u[,1])
plot(wavelengths, -svd.temp$u[,2])
plot(wavelengths, -svd.temp$u[,3])

plot(temperatures, svd.temp$v[,1]*svd.temp$d[1])
plot(temperatures, svd.temp$v[,2])
plot(temperatures, svd.temp$v[,3])

y1 <- svd.temp$v[,1]*svd.temp$d[1]
y2 <- svd.temp$v[,2]*svd.temp$d[2]
y3 <- svd.temp$v[,3]*svd.temp$d[3]
y4 <- svd.temp$v[,4]*svd.temp$d[4]
y5 <- svd.temp$v[,5]*svd.temp$d[5]
y6 <- svd.temp$v[,6]*svd.temp$d[6]
x <- temperatures

plot(x, -y1)
plot(x, y2)
plot(x, y3)
plot(x, y4)
plot(x, y5)
plot(x, y6)

plot(svd.temp$v%*%diag(svd.temp$d)[,1])

svd.temp$u%*%diag(svd.temp$d)%*%t(svd.temp$v)

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

#############################
### Minimization with 1, 2, 3, 4, or 5 vectors
############

one.state.1 <- function(z)
{
  Sn1<-z[1]
  Su1<-z[2]
  dH<-z[3]
  Tm<-z[4]
  sum((y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2)  ### How can we do this in a reduced way.....
  ### One response vector that is k*length(x).  Then 
}

r1 <- optim(c(0,0,-20000, 350), one.state.1, method='L-BFGS-B')
r1
plot(x,y1)
lines(x, SingleStateMeltingFunction(x, r1$par[1], r1$par[2], r1$par[3], r1$par[4]))

one.state.2 <- function(z)
{
  Sn1<-z[1]
  Su1<-z[2]
  Sn2<-z[3]
  Su2<-z[4]
  dH<-z[5]
  Tm<-z[6]
  sum((y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2) +
    sum((y2 - SingleStateMeltingFunction(x, Sn2, Su2, dH, Tm))^2)
}

r2 <- optim(c(0,0,0,0,-25000, 350), one.state.2, method='L-BFGS-B')
r2
plot(x,y1)
lines(x, SingleStateMeltingFunction(x, r2$par[1], r2$par[2], r2$par[5], r2$par[6]))
plot(x,y2)
lines(x, SingleStateMeltingFunction(x, r2$par[3], r2$par[4], r2$par[5], r2$par[6]))

### A fit of one-transition to 3, 4, and 5 vectors

one.state.3 <- function(z)
{
  Sn1 <- z[1]
  Su1 <- z[2]
  Sn2 <- z[3]
  Su2 <- z[4]
  Sn3 <- z[5]
  Su3 <- z[6]
  dH <- z[7]
  Tm <- z[8]
  sum((y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2) +
    sum((y2 - SingleStateMeltingFunction(x, Sn2, Su2, dH, Tm))^2) +
    sum((y3 - SingleStateMeltingFunction(x, Sn3, Su3, dH, Tm))^2)
}

r2 <- optim(c(0,0,0,0,0,0,-25000, 350), one.state.3, method='L-BFGS-B')
r2
plot(x,y1)
lines(x, SingleStateMeltingFunction(x, r2$par[1], r2$par[2], r2$par[7], r2$par[8]))
plot(x,y2)
lines(x, SingleStateMeltingFunction(x, r2$par[3], r2$par[4], r2$par[7], r2$par[8]))
plot(x,y3)
lines(x, SingleStateMeltingFunction(x, r2$par[5], r2$par[6], r2$par[7], r2$par[8]))

##############################
### Two State Functions ###
##############################

TwoStateMeltingFunction <- function(T,Sn,S1,Su,dH1,Tm1,dH2,Tm2) #Native to Intermediate 1 to Unfolded
{
  (Sn + S1*exp(-((dH1/R)*((1/Tm1) - (1/T)))) + 
     Su*exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T)))) +
             exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))
}

#############################
### Minimization with 1, 2, 3, 4, or 5 vectors
############

two.state.1 <- function(z)
{
  Sn1<-z[1]
  S11<-z[2]
  Su1<-z[3]
  dH1<-z[4]
  Tm1<-z[5]
  dH2<-z[6]
  Tm2<-z[7]
  sum((y1 - TwoStateMeltingFunction(x, Sn1, S11, Su1, dH1, Tm1, dH2, Tm2))^2)  ### How can we do this in a reduced way.....
  ### One response vector that is k*length(x).  Then 
}

r1 <- optim(c(0,0,0,-25000, 310, -30000, 350), two.state.1, method='L-BFGS-B')
r1
plot(x,y1)
lines(x, TwoStateMeltingFunction(x, r1$par[1], r1$par[2], r1$par[3], r1$par[4], r1$par[5], r1$par[6], r1$par[7]))

two.state.2 <- function(z)
{
  Sn1<-z[1]
  S11<-z[2]
  Su1<-z[3]
  Sn2<-z[4]
  S12<-z[5]
  Su2<-z[6]
  dH1<-z[7]
  Tm1<-z[8]
  dH2<-z[9]
  Tm2<-z[10]
  sum((y1 - TwoStateMeltingFunction(x, Sn1, S11, Su1, dH1, Tm1, dH2, Tm2))^2) +
    sum((y2 - TwoStateMeltingFunction(x, Sn2, S12, Su2, dH1, Tm1, dH2, Tm2))^2)
}

r2 <- optim(c(0,0,0,0,0,0, -24000, 323, -31000, 341), two.state.2, method='L-BFGS-B')
r2
r3 <- optim(r2$par, two.state.2, method='L-BFGS-B')
r3

for(j in 1:1000)
  r3 <- optim(r3$par, two.state.2, method='L-BFGS-B')

r3
r2<-r3
plot(x,y1)
lines(x, TwoStateMeltingFunction(x, r2$par[1], r2$par[2], r2$par[3], r2$par[7], r2$par[8], r2$par[9], r2$par[10]))
plot(x,y2)
lines(x, TwoStateMeltingFunction(x, r2$par[4], r2$par[5], r2$par[6], r2$par[7], r2$par[8], r2$par[9], r2$par[10]))


#### this is a fit of 3-vectors to a two-transition model

two.state.3 <- function(z)
{
  Sn1<-z[1]
  S11<-z[2]
  Su1<-z[3]
  Sn2<-z[4]
  S12<-z[5]
  Su2<-z[6]
  Sn3<-z[7]
  S13<-z[8]
  Su3<-z[9]
  dH1<-z[10]
  Tm1<-z[11]
  dH2<-z[12]
  Tm2<-z[13]
  sum((y1 - TwoStateMeltingFunction(x, Sn1, S11, Su1, dH1, Tm1, dH2, Tm2))^2) +
    sum((y2 - TwoStateMeltingFunction(x, Sn2, S12, Su2, dH1, Tm1, dH2, Tm2))^2) +
    sum((y3 - TwoStateMeltingFunction(x, Sn3, S13, Su3, dH1, Tm1, dH2, Tm2))^2)
}

r2 <- optim(c(0,0,0,0,0,0,0,0,0, -24000, 323, -31000, 341), two.state.3, method='L-BFGS-B')
r2
r3 <- optim(r2$par, two.state.3, method='L-BFGS-B')
r3

for(j in 1:1000)
  r3 <- optim(r3$par, two.state.3, method='L-BFGS-B')

r3
r2<-r3
plot(x,y1)
lines(x, TwoStateMeltingFunction(x, r2$par[1], r2$par[2], r2$par[3], r2$par[10], r2$par[11], r2$par[12], r2$par[13]))
plot(x,y2)
lines(x, TwoStateMeltingFunction(x, r2$par[4], r2$par[5], r2$par[6], r2$par[10], r2$par[11], r2$par[12], r2$par[13]))
plot(x,y3)
lines(x, TwoStateMeltingFunction(x, r2$par[7], r2$par[8], r2$par[9], r2$par[10], r2$par[11], r2$par[12], r2$par[13]))
