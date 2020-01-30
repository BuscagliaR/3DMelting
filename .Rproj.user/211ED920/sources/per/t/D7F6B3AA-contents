First.Matrix <- read.table('22wt_Melt_2_CD.txt')
Temperatures = seq(10, 102, 2) ### GET THIS FROM PUB
Wavelengths = seq(300, 600, 2) ### GET THIS FROM PUB
plot(First.Matrix[,1])
plot(First.Matrix[,30])

plot(as.numeric(First.Matrix[70,]))


SVD.temp <- svd(First.Matrix)

plot(SVD.temp$d)

plot(SVD.temp$u[,1])

plot(Temperatures, SVD.temp$v[,1])


SingleStateMeltingFunction <- function(T, Sn, Su, dH, Tm)
{
  R = 1.986 ### Gas constant
  (Sn + Su*exp(-(dH/R)*(1/Tm - 1/T)))/(1+exp(-(dH/R)*(1/Tm - 1/T)))
}

formula(SingleStateMeltingFunction(T, 'Sn1', 'Su1', 'dH', 'Tm'))

### formula y ~ x1 + x2 + x3
### formula y ~ 


#### Simulation from https://www.r-bloggers.com/first-steps-with-non-linear-regression-in-r/

#simulate some data
set.seed(20160227)
x<-seq(0,50,1)
y<-((runif(1,10,20)*x)/(runif(1,0,10)+x))+rnorm(51,0,1)
#for simple models nls find good starting values for the parameters even if it throw a warning
m<-nls(y~a*x/(b+x))

summary(m)

#get some estimation of goodness of fit
cor(y,predict(m))

plot(x,y)

R = 1.986
y = SVD.temp$v[,1]
x = Temperatures+273.15

plot(x,y)

m<-nls(y~Sn + Su*exp(-(dH/R)*(1/Tm - 1/x))/(1+exp(-(dH/R)*(1/Tm - 1/x))),
       start=list(Sn=-0.2, Su=0.1, dH = 2000, Tm = 350))
summary(m)

(Sn + Su*exp(-(dH/R)*(1/Tm - 1/x)))/(1+exp(-(dH/R)*(1/Tm - 1/x)))
