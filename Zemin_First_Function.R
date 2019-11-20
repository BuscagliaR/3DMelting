#Set the constant R value
R = 1.986

SingleStateMeltingFunction <- function(T, Sn, Su, dH, Tm)
{
  (Sn + Su*exp(-(dH/R)*(1/Tm - 1/T)))/(1+exp(-(dH/R)*(1/Tm - 1/T)))
}

T.grid = seq(0, 100, 2) #Temperatures from 20 to 100 by 2 degrees
y1 <- SingleStateMeltingFunction(T = T.grid, Sn = 6, Su = 15, dH = -200, Tm = 20)
y2 <- SingleStateMeltingFunction(T = T.grid, Sn = 15, Su = 0, dH = -2000, Tm = 50)

plot(T.grid, y1)
plot(T.grid, y2)
plot(T.grid, y1 + y2)


SingleStateFractionSn <- function(T,dH,Tm)
{
  1/(1 + exp(-((dH/R)*((1/Tm) - (1/T)))))
}

T.grid = seq(0, 100, 2) #Temperatures from 20 to 100 by 2 degrees
y1 <- SingleStateFractionSn(T = T.grid, dH = -200, Tm = 20)
y2 <- SingleStateFractionSn(T = T.grid, dH = -2000, Tm = 50)

plot(T.grid, y1)
plot(T.grid, y2)
plot(T.grid, y1 + y2)

SingleleStateFractionSu <- function(T,dH,Tm)
{
  (exp(-((dH/R)*((1/Tm) - (1/T)))))/(1 + exp(-((dH/R)*((1/Tm) - (1/T)))))
}

y3 <- SingleleStateFractionSu(T = T.grid, dH = -2000, Tm = 50)
plot(T.grid,y2)
points(T.grid,y3)



TwoStateMeltingFunction <- function(T,Sn,S1,Su,dH1,Tm1,dH2,Tm2) #Native to Intermediate 1 to Unfolded
{
  (Sn + S1*exp(-((dH1/R)*((1/Tm1) - (1/T)))) + 
      Su*exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T)))) +
                                                                                 exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))
  
}

y4 <- TwoStateMeltingFunction(T = T.grid,Sn = 100,S1 = 25,Su = 75, dH1 = -2000, Tm1 = 40, dH2 = -3000, Tm2 = 65)
plot(T.grid,y4)


TwoStateFractionSn <- function(T,dH1,Tm1,dH2,Tm2)
{
  1/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T))))) + exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T))))
}

TwoStateFractionS1 <- function(T,dH1,Tm1,dH2,Tm2)
{
  exp(-((dH1/R)*((1/Tm1) - (1/T))))/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T)))) + exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))
}

TwoStateFractionSu <- function(T,dH1,Tm1,dH2,Tm2)
{
  (exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T)))) + exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))))
  
}

K1 = function(T,dH1,Tm1)
{
  (exp((-dH1/R)*(1/Tm1 - 1/T)))
}

T.grid <- seq(20, 100, 1)
plot(T.grid, K1(T.grid, -2000, 45))

K2 = function(T,dH2,Tm2)
{
  (exp((-dH2/R)*(1/Tm2 - 1/T)))
}

K0 <- function(T,dH1,Tm1,dH2,Tm2)
{
  (K1(T,dH1,Tm1)/K2(T,dH2,Tm2))
}

plot(T.grid, K0(T.grid, -2000, 40, -2500, 60))

Denom <- 1+K0(T,dH1,Tm1,dH2,Tm2) + K0(T,dH1,Tm1,dH2,Tm2)*K2(T,dH2,Tm2)

ParallelFractionSn1 <- function(T,dH1,Tm1,dH2,Tm2)
{
  (1/Denom)
}

ParallelFractionSn2 <- function(T,dH1,Tm1,dH2,Tm2)
{
  (K0(T,dH1,Tm1,dH2,Tm2)/Denom)
}

ParallelFractionSu <- function(T,dH1,Tm1,dH2,Tm2)
{
  (K1(T,dH1,Tm1)/Denom)
}

ParallelSpecies <- function(T,Sn1,Sn2,Su,dH1,Tm1,dH2,Tm2)
{
  (Sn1*ParallelFractionSn1 + Sn2*ParallelFractionSn2+ Su*ParallelFractionSu)
}

ThreeStateMeltingFunction <- function(T,Sn,S1,S2,Su,dH1,Tm1,dH2,Tm2,dH3,Tm3)
{
(Sn + S1*exp(-((dH1/R)*((1/Tm1) - (1/T))))) + S2*exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T)))) + Su*exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T))))*exp(-((dH3/R)*((1/Tm3) - (1/T))))/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T)))) + exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T))))) + exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T))))*exp(-((dH3/R)*((1/Tm3) - (1/T))))

}

ThreestateFractionSn <- function(T,dH1,Tm1,dH2,Tm2,dH3,Tm3)
{
  
1/(1 + exp(-((dH1/R)*((1/Tm1) - (1/T)))) + 
      exp(-((dH1/R)*((1/Tm1) - (1/T))))*
      exp(-((dH2/R)*((1/Tm2) - (1/T))))+ 
      exp(-((dH1/R)*((1/Tm1) - (1/T))))*
      exp(-((dH2/R)*((1/Tm2) - (1/T))))*
      exp(-((dH3/R)*((1/Tm3) - (1/T)))))
}

ThreestateFractionS1<- function(T,dH1,Tm1,dH2,Tm2,dH3,Tm3)
{
  exp(-((dH1/R)*((1/Tm1) - (1/T))))/(1 + 
                                       exp(-((dH1/R)*((1/Tm1) - (1/T)))) + 
                                       exp(-((dH1/R)*((1/Tm1) - (1/T))))*
                                       exp(-((dH2/R)*((1/Tm2) - (1/T))))+ 
                                       exp(-((dH1/R)*((1/Tm1) - (1/T))))*
                                       exp(-((dH2/R)*((1/Tm2) - (1/T))))*
                                       exp(-((dH3/R)*((1/Tm3) - (1/T)))));
  
}

ThreestateFractionS2<- function(T,dH1,Tm1,dH2,Tm2,dH3,Tm3)
{(exp(-((dH1/R)*((1/Tm1) - (1/T))))*
    exp(-((dH2/R)*((1/Tm2) - (1/T)))))/(1 + 
                                          exp(-((dH1/R)*((1/Tm1) - (1/T)))) + 
                                          exp(-((dH1/R)*((1/Tm1) - (1/T))))*
                                          exp(-((dH2/R)*((1/Tm2) - (1/T)))) + 
                                          exp(-((dH1/R)*((1/Tm1) - (1/T))))*
                                          exp(-((dH2/R)*((1/Tm2) - (1/T))))*
                                          exp(-((dH3/R)*((1/Tm3) - (1/T)))));
}

ThreestateFractionSu <- function(T,dH1,Tm1,dH2,Tm2,dH3,Tm3)
{
  exp(-((dH1/R)*((1/Tm1) - (1/T))))*exp(-((dH2/R)*((1/Tm2) - (1/T))))*
    exp(-((dH3/R)*((1/Tm3) - (1/T)))) / (1 + 
                                          exp(-((dH1/R)*((1/Tm1) - (1/T)))) + 
                                          exp(-((dH1/R)*((1/Tm1) - (1/T))))*
                                          exp(-((dH2/R)*((1/Tm2) - (1/T)))) + 
                                          exp(-((dH1/R)*((1/Tm1) - (1/T))))*
                                          exp(-((dH2/R)*((1/Tm2) - (1/T))))*
                                          exp(-((dH3/R)*((1/Tm3) - (1/T)))))
}

K1Int <- function(T, dH1, Tm1)
{
  exp(-(dH1/R)*(1/Tm1 - 1/T))
}

K2Int <- function(T, dH2,Tm2)
{
  exp(-(dH2/R)*(1/Tm2 - 1/T))
}

KUInt <- function(T,dH3,Tm3)
{ 
  exp(-(dH3/R)*(1/Tm3 - 1/T))
}

DenomInt <- function(T, dH1, Tm1, dH2, Tm2, dH3, Tm3)
{
  1 + K1Int(T, dH1, Tm1) / K2Int(T, dH2, Tm2) + K1Int(T, dH1, Tm1) + K1Int(T, dH1, Tm1)*KUInt(T, dH3, Tm3)
}

T.grid<-seq(20, 100, 1)
plot(T.grid, y=DenomInt(T.grid, -2000, 40, -2500, 50, -3000, 60))

#put in parameters in each function 
ParallelIntFractionSn1 <- function(T, dH1, Tm1, dH2, Tm2, dH3, Tm3)
{
  1 / DenomInt(T, dH1, Tm1, dH2, Tm2, dH3, Tm3)
}

ParallelIntFractionSn2 <- function(T, dH1, Tm1, dH2, Tm2, dH3, Tm3)
{
  (K1Int(T,dH1,Tm1) / K2Int(T,dH2,Tm2) / DenomInt(T,dH1,Tm1,dH2,Tm2,dH3,Tm3))
}

ParallelIntFractionSI <- function(T, dH1,Tm1, dH2, Tm2, dH3, Tm3)
{
  K1Int(T,dH1,Tm1) / DenomInt(T,dH1,Tm1,dH2,Tm2,dH3,Tm3)
}

ParallelIntFractionSU <- function(T, dH1, Tm1, dH2, Tm2, dH3, Tm3)
{
  KUInt(T,dH3,Tm3) * ParallelIntFractionSI(T, dH1,Tm1, dH2, Tm2, dH3, Tm3)
}

ParallelIntermediateFunction <- function(T, Sn1, Sn2, S1, Su, dH1, Tm1, dH2, Tm2, dH3, Tm3)
{
  Sn1 * ParallelIntFractionSn1(T, dH1, Tm1, dH2, Tm2, dH3, Tm3) + Sn2 * ParallelIntFractionSn2(T, dH1, Tm1, dH2, Tm2, dH3, Tm3) +S1 * ParallelIntFractionSI(T, dH1,Tm1, dH2, Tm2, dH3, Tm3) + Su * ParallelIntFractionSU(T, dH1, Tm1, dH2, Tm2, dH3, Tm3)
}
T.grid <- seq(20, 100, 1)
y <- ParallelIntermediateFunction(T.grid, 500, 800, 200, 0, -2000, 40, -2500, 50, -3000, 70)
plot(T.grid, y)
''