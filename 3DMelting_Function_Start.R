### As discussed, it seems like R may require these to be in a slightly different form
### While we are working on that, I wanted to confirm the functions were working.


### Setup ####

### Data

raw.data <- read.table('2GKU_Melt_2_CD.txt')
temperatures <- seq(4, 94, 2) + 273.15
wavelengths <- seq(220, 340, 1)

### Visualize, maybe a 3d plot?

plot(wavelengths, raw.data[,1])
plot(temperatures, raw.data[15,])

### Need to run 'estimated' parameters simulation as in Mathematica

### 
svd.temp <- svd(raw.data)

plot(wavelengths, -svd.temp$u[,1])
plot(wavelengths, -svd.temp$u[,2])
plot(wavelengths, -svd.temp$u[,3])

plot(temperatures, svd.temp$v[,1]*svd.temp$d[1])
plot(temperatures, svd.temp$v[,2])
plot(temperatures, svd.temp$v[,3])

y1 <- svd.temp$v[,1]*svd.temp$d[1]
y2 <- svd.temp$v[,2]*svd.temp$d[2]
y3 <- svd.temp$v[,3]*svd.temp$d[1]
y4 <- svd.temp$v[,4]*svd.temp$d[2]
y5 <- svd.temp$v[,5]*svd.temp$d[5]
x <- temperatures

plot(svd.temp$v%*%diag(svd.temp$d)[,1])

svd.temp$u%*%svd.temp$d%*%t(svd.temp$v)

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
r1$value
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

r2 <- optim(c(0,0,0,0,0,0, -24000, 341, -31000, 323), two.state.2, method='L-BFGS-B')
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
