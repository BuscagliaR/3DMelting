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

y1 <- svd.temp$v[,1]*svd.temp$d[1]
y2 <- svd.temp$v[,2]*svd.temp$d[2]
x <- temperatures

fun.test.1 <- function(Sn1, Su1, dH, Tm)
{
  (y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2
}

optim(c(1,1,-20000, 350), fun.test.1)

fun.test.3 <- function(z)
{
  Sn1<-z[1]
  Su1<-z[2]
  dH<-z[3]
  Tm<-z[4]
  sum((y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2)
}

r1 <- optim(c(0,0,-15000, 290), fun.test.3, method='L-BFGS-B')
# r1 <- optim(c(0,0,-20000, 350), fun.test.3)
r1

r1$par

plot(x,y1)
lines(x, SingleStateMeltingFunction(x, r1$par[1], r1$par[2], r1$par[3], r1$par[4]))




### optimize over two vectors with different signals

fun.test.4 <- function(z)
{
  Sn1<-z[1]
  Su1<-z[2]
  Sn2<-z[3]
  Su2<-z[4]
  dH<-z[5]
  Tm<-z[6]
  out <- sum((y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2) +
    sum((y2 - SingleStateMeltingFunction(x, Sn2, Su2, dH, Tm))^2)
}

r2 <- optim(c(0,0,0,0,-15000, 290), fun.test.4, method='L-BFGS-B')
# r1 <- optim(c(0,0,-20000, 350), fun.test.3)
r2

r1$par

plot(x,y1)
lines(x, SingleStateMeltingFunction(x, r2$par[1], r2$par[2], r2$par[5], r2$par[6]))

plot(x,y2)
lines(x, SingleStateMeltingFunction(x, r2$par[3], r2$par[4], r2$par[5], r2$par[6]))

### Woo its working!


fun.test.4 <- function(z)
{
  Sn1<-z[1]
  Su1<-z[2]
  Sn2<-z[3]
  Su2<-z[4]
  dH<-z[5]
  Tm<-z[6]
  sum((y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2)
}


fun.test.2 <- function(Sn1, Su1, Sn2, Su2, dH, Tm)
{
  (y1 - SingleStateMeltingFunction(x, Sn1, Su1, dH, Tm))^2 + (y2 - SingleStateMeltingFunction(x, Sn2, Su2, dH, Tm))^2
}

optim(c(1,1,1,1,-20000, 350), fun.test.2)

optimize(fun.test)

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
  
  
  
  
  
  
  
  
  
  
  

#### HELP 


## Combinatorial optimization: Traveling salesman problem
library(stats) # normally loaded

eurodistmat <- as.matrix(eurodist)

distance <- function(sq) {  # Target function
  sq2 <- embed(sq, 2)
  sum(eurodistmat[cbind(sq2[,2], sq2[,1])])
}

genseq <- function(sq) {  # Generate new candidate sequence
  idx <- seq(2, NROW(eurodistmat)-1)
  changepoints <- sample(idx, size = 2, replace = FALSE)
  tmp <- sq[changepoints[1]]
  sq[changepoints[1]] <- sq[changepoints[2]]
  sq[changepoints[2]] <- tmp
  sq
}

sq <- c(1:nrow(eurodistmat), 1)  # Initial sequence: alphabetic
distance(sq)
# rotate for conventional orientation
loc <- -cmdscale(eurodist, add = TRUE)$points
x <- loc[,1]; y <- loc[,2]
s <- seq_len(nrow(eurodistmat))
tspinit <- loc[sq,]

plot(x, y, type = "n", asp = 1, xlab = "", ylab = "",
     main = "initial solution of traveling salesman problem", axes = FALSE)
arrows(tspinit[s,1], tspinit[s,2], tspinit[s+1,1], tspinit[s+1,2],
       angle = 10, col = "green")
text(x, y, labels(eurodist), cex = 0.8)

set.seed(123) # chosen to get a good soln relatively quickly
res <- optim(sq, distance, genseq, method = "SANN",
             control = list(maxit = 30000, temp = 2000, trace = TRUE,
                            REPORT = 500))
res  # Near optimum distance around 12842

tspres <- loc[res$par,]
plot(x, y, type = "n", asp = 1, xlab = "", ylab = "",
     main = "optim() 'solving' traveling salesman problem", axes = FALSE)
arrows(tspres[s,1], tspres[s,2], tspres[s+1,1], tspres[s+1,2],
       angle = 10, col = "red")
text(x, y, labels(eurodist), cex = 0.8)
