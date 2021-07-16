################################################################################
#R code "Divergent pathways of nitrogen-fixing trees through succession ######## 
#depend on starting nitrogen supply and priority effects" #######################
################################################################################

################################################################################
#Base Parameter Values (Table S1) ##############################################
################################################################################

I <- 4 #Abiotic N input flux [kgN/ha/y]
m <- 0.01 #Net N mineralization rate [1/y]
k <- 5 #Available N loss rate [1/y]
uF <- 0.2 #Fixer turnover rate [1/y] 
u0 <- 0.2 #Nonfixer turnover rate [1/y]
wF <- 50 #Fixer N use efficiency [kgC/kgN]
w0 <- 50 #Nonfixer N use efficiency [kgC/kgN]
vF <- 0.05 #Fixer N uptake parameter [ha/kgC/y] 
v0 <- 0.07 #Nonfixer N uptake parameter [ha/kgC/y]
f <- 0.01 #N-fixation parameter [kgN/kgC/y]
phi <- 0.001 #Unavailable N loss rate [1/y]
bF <- 2.1 #Fixer non-N-limited growth parameter 1 [1/y]
b0 <- 2 #Nonfixer non-N-limited growth parameter 1 [1/y]
gF <- 0.00025 #Fixer non-N-limited growth parameter 2 [ha/kgC]
g0 <- 0.0002 #Nonfixer non-N-limited growth parameter 2 [ha/kgC]

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#Set up time
t0 <- 1
tf <- 10000
tstep <- 1
ts <- seq(t0,tf,tstep)

#Install package to numerically solve our system of equations
install.packages("deSolve")
library(deSolve)

################################################################################
#Function definitions ##########################################################
################################################################################

#timescale approximation model (BF, B0, A; facultative SNF)
static.D.fac <- function(t,x,parms){
  with(as.list(c(parms,x)),{
    Ff <- min(max(bF/(wF*(1 + gF*(x[1] + x[2]))) - vF*x[3],0),f) #Equation 7
    gf <- min(wF*vF*x[3] + wF*Ff,bF/(1 + gF*(x[1] + x[2]))) #Equation 5
    g0 <- min(w0*v0*x[3],b0/(1 + g0*(x[1] + x[2]))) #Equation 6
    dBFdt <- x[1]*(gf - uF) #Equation 1
    dB0dt <- x[2]*(g0 - u0) #Equation 2
    dAdt <- I + m*D - k*x[3] - x[1]/wF*(gf - wF*Ff) - x[2]/w0*g0 #Equation 4
    return(list(c(dBFdt,dB0dt,dAdt)))
  })
}

#timescale approximation model (BF, B0, A; obligate SNF)
static.D.obl <- function(t,x,parms){
  with(as.list(c(parms,x)),{
    gf <- min(wF*vF*x[3] + wF*f,bF/(1 + gF*(x[1] + x[2]))) #Equation 5
    g0 <- min(w0*v0*x[3],b0/(1 + g0*(x[1] + x[2]))) #Equation 6
    dBFdt <- x[1]*(gf - uF) #Equation 1
    dB0dt <- x[2]*(g0 - u0) #Equation 2
    dAdt <- I + m*D - k*x[3] - x[1]/wF*(gf - wF*f) - x[2]/w0*g0 #Equation 4
    return(list(c(dBFdt,dB0dt,dAdt)))
  })
}

#full model (facultative SNF)
dynamic.D.fac <- function(t,x,parms){
  with(as.list(c(parms,x)),{
    Ff <- min(max(bF/(wF*(1 + gF*(x[1] + x[2]))) - vF*x[4],0),f) #Equation 7
    gf <- min(wF*vF*x[4] + wF*Ff,bF/(1 + gF*(x[1] + x[2]))) #Equation 5
    g0 <- min(w0*v0*x[4],b0/(1 + g0*(x[1] + x[2]))) #Equation 6
    dBFdt <- x[1]*(gf - uF) #Equation 1
    dB0dt <- x[2]*(g0 - u0) #Equation 2
    dDdt <- uF*x[1]/wF + u0*x[2]/w0 - x[3]*(m + phi) #Equation 3
    dAdt <- I + m*x[3] - k*x[4] - x[1]/wF*(gf - wF*Ff) - x[2]/w0*g0 #Equation 4
    return(list(c(dBFdt,dB0dt,dDdt,dAdt)))
  })
}

#full model (obligate SNF)
dynamic.D.obl <- function(t,x,parms){
  with(as.list(c(parms,x)),{
    gf <- min(wF*vF*x[4] + wF*f,bF/(1 + gF*(x[1] + x[2]))) #Equation 5
    g0 <- min(w0*v0*x[4],b0/(1 + g0*(x[1] + x[2]))) #Equation 6
    dBFdt <- x[1]*(gf - uF) #Equation 1
    dB0dt <- x[2]*(g0 - u0) #Equation 2
    dDdt <- uF*x[1]/wF + u0*x[2]/w0 - x[3]*(m + phi) #Equation 3
    dAdt <- I + m*x[3] - k*x[4] - x[1]/wF*(gf - wF*f) - x[2]/w0*g0 #Equation 4
    return(list(c(dBFdt,dB0dt,dDdt,dAdt)))
  })
}

################################################################################
#Simulations in figure 3 #######################################################
################################################################################

#solve for critical S values for figure 3A,B,E,F,G,H (SF<Sb<Sa<SF.star<S0)
print(Sa <- (k + v0*(bF - uF)/(gF*uF))*u0/(w0*v0)) #Equation 11
print(Sb <- min((k + vF*(bF - uF)/(gF*uF))*u0/(w0*v0),k*u0/(w0*v0)+uF*((bF - uF)/(gF*uF))/wF)) #Equation 12
print(SF <- (k + vF*(bF - uF)/(gF*uF))*(uF - wF*f)/(wF*vF)) #Equation 13
print(S0 <- (k + v0*(b0 - u0)/(g0*u0))*u0/(w0*v0)) #Equation 14
print(SF.star <- (k + vF*(bF - uF)/(gF*uF))*(uF)/(wF*vF)) #Equation 13 with downregulated SNF (facultative SNF)

#find corresponding D values to critical S values for figure 3A,B,E,F,G,H
#D <- (S - I)/m
print(Da <- (Sa - I)/m)
print(Db <- (Sb - I)/m)
print(DF <- (SF - I)/m)
print(D0 <- (S0 - I)/m)
print(DF.star <- (SF.star - I)/m)

#equilibrium A values under N-limitation (equation 8)
AbarBF <- (uF - wF*f)/(vF*wF)
AbarB0 <- u0/(v0*w0) 
AbarBF;AbarB0

#equilibrium plant biomass values under R2-limitation (equation 9)
BFbar <- (bF - uF)/(gF*uF)
B0bar <- (b0 - u0)/(g0*u0)
BFbar;B0bar

#Competitive Exclusion of Nonfixers (figure 3A,B)
#simulate successional trajectories
D <- 1000 #S<Sb
sim.3A.static <- lsoda(c(8000,2000,0),ts,static.D.fac,parms)
sim.3A.dynamic <- lsoda(c(8000,2000,D,0),ts,dynamic.D.fac,parms)

sim.3B.static <- lsoda(c(2000,8000,AbarB0),ts,static.D.fac,parms)
sim.3B.dynamic <- lsoda(c(2000,8000,D,AbarB0),ts,dynamic.D.fac,parms)

#Figure 3A
par(pty="s")
plot(sim.3A.static[,1],sim.3A.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3A.static[,1],sim.3A.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3A.dynamic[,1],sim.3A.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3A.dynamic[,1],sim.3A.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="N-fixers exclude",line=1.5,cex.main=1.5)
title(main="nonfixers",line=0.2,cex.main=1.5)

#Figure 3B
plot(sim.3B.static[,1],sim.3B.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3B.static[,1],sim.3B.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3B.dynamic[,1],sim.3B.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3B.dynamic[,1],sim.3B.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="N-fixers exclude",line=1.5,cex.main=1.5)
title(main="nonfixers",line=0.2,cex.main=1.5)

#Coexistence (figure 3E,F)
#simulate successional trajectories
D <- 14000 #Sb<S<Sa
sim.3E.static <- lsoda(c(8000,2000,0),ts,static.D.fac,parms)
sim.3E.dynamic <- lsoda(c(8000,2000,D,0),ts,dynamic.D.fac,parms)

sim.3F.static <- lsoda(c(2000,8000,AbarB0),ts,static.D.fac,parms)
sim.3F.dynamic <- lsoda(c(2000,8000,D,AbarB0),ts,dynamic.D.fac,parms)

#Figure 3E
plot(sim.3E.static[,1],sim.3E.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3E.static[,1],sim.3E.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3E.dynamic[,1],sim.3E.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3E.dynamic[,1],sim.3E.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="Coexistence",line=1,cex.main=1.5)

#Figure 3F
plot(sim.3F.static[,1],sim.3F.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3F.static[,1],sim.3F.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3F.dynamic[,1],sim.3F.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3F.dynamic[,1],sim.3F.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="Coexistence",line=1,cex.main=1.5)

#Competitive Exclusion of N-fixers (figure 3G,H)
#simulate successional trajectories
D <- 30000 #Sa<S
sim.3G.static <- lsoda(c(8000,2000,0),ts,static.D.fac,parms)
sim.3G.dynamic <- lsoda(c(8000,2000,D,0),ts,dynamic.D.fac,parms)

sim.3H.static <- lsoda(c(2000,8000,AbarB0),ts,static.D.fac,parms)
sim.3H.dynamic <- lsoda(c(2000,8000,D,AbarB0),ts,dynamic.D.fac,parms)

#Figure 3G
plot(sim.3G.static[,1],sim.3G.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3G.static[,1],sim.3G.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3G.dynamic[,1],sim.3G.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3G.dynamic[,1],sim.3G.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="Nonfixers exclude",line=1.5,cex.main=1.5)
title(main="N-fixers",line=0.2,cex.main=1.5)

#Figure 3H
plot(sim.3H.static[,1],sim.3H.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3H.static[,1],sim.3H.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3H.dynamic[,1],sim.3H.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3H.dynamic[,1],sim.3H.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="Nonfixers exclude",line=1.5,cex.main=1.5)
title(main="N-fixers",line=0.2,cex.main=1.5)

#Priority effects (figure 3C,D)
#Parameter change
I <- 8
wF <- 40
vF <- 0.1
phi <- 0.0005

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#solve for critical S values for figure 3C,D (SF<Sa<S0<SF.star<Sb;Sa & S0 remain unchanged)
print(Sb <- min((k + vF*(bF - uF)/(gF*uF))*u0/(w0*v0),k*u0/(w0*v0)+uF*((bF - uF)/(gF*uF))/wF)) #Equation 12
print(SF <- (k + vF*(bF - uF)/(gF*uF))*(uF - wF*f)/(wF*vF)) #Equation 13
print(SF.star <- (k + vF*(bF - uF)/(gF*uF))*(uF)/(wF*vF)) #Equation 13 with downregulated SNF (facultative SNF)

#find corresponding D values to critical S values for figure 3C,D
#D <- (S - I)/m
print(Db <- (Sb - I)/m)
print(DF <- (SF - I)/m)
print(DF.star <- (SF.star - I)/m)

#equilibrium A values under N-limitation (equation 8; AbarB0 does not change with these parameters)
AbarBF <- (uF - wF*f)/(vF*wF)
AbarBF;AbarB0

#equilibrium plant biomass values under R2-limitation remain unchanged (equation 9)

#simulate successional trajectories
D <- 16000 #Sa<S<Sb
sim.3C.static <- lsoda(c(8000,2000,0),ts,static.D.fac,parms)
sim.3C.dynamic <- lsoda(c(8000,2000,D,0),ts,dynamic.D.fac,parms)

sim.3D.static <- lsoda(c(2000,8000,AbarB0),ts,static.D.fac,parms)
sim.3D.dynamic <- lsoda(c(2000,8000,D,AbarB0),ts,dynamic.D.fac,parms)

#Figure 3C
plot(sim.3C.static[,1],sim.3C.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3C.static[,1],sim.3C.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3C.dynamic[,1],sim.3C.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3C.dynamic[,1],sim.3C.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="Priority effects",line=1,cex.main=1.5)

#Figure 3D
plot(sim.3D.static[,1],sim.3D.static[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,xlim=c(0,600),ylim=c(0,60000),las=1)
lines(sim.3D.static[,1],sim.3D.static[,3],col="blue",lty=2,lwd=2)
lines(sim.3D.dynamic[,1],sim.3D.dynamic[,2],col="red",lty=1,lwd=0.5)
lines(sim.3D.dynamic[,1],sim.3D.dynamic[,3],col="blue",lty=2,lwd=0.5)
mtext(expression('Biomass (kg C ha'^-1*')'),side=2,line=4,cex=1.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)
title(main="Priority effects",line=1,cex.main=1.5)

################################################################################
#Simulations in figure S2 ######################################################
################################################################################

#Set parameters as in Table S1 except for wF,w0, & uF which are specific to this figure
I <- 4
uF <- 0.25
wF <- 25
w0 <- 60
vF <- 0.05
phi <- 0.001

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#solve for critical S values for figure S2 (SF<Sb<Sa<S0)
print(Sa <- (k + v0*(bF - uF)/(gF*uF))*u0/(w0*v0)) #Equation 11
print(Sb <- k*u0/(w0*v0)+(uF-wF*f)*((bF - uF)/(gF*uF))/wF) #Equation 12
print(SF <- (k + vF*(bF - uF)/(gF*uF))*(uF - wF*f)/(wF*vF)) #Equation 13
print(S0 <- (k + v0*(b0 - u0)/(g0*u0))*u0/(w0*v0)) #Equation 14

#find corresponding D values to critical S values for figure S2
#D <- (S - I)/m
print(Da <- (Sa - I)/m)
print(Db <- (Sb - I)/m)
print(DF <- (SF - I)/m)
print(D0 <- (S0 - I)/m)

#equilibrium A values under N-limitation (equation 8)
AbarBF <- (uF - wF*f)/(vF*wF)
AbarB0 <- u0/(v0*w0) 
AbarBF;AbarB0

#equilibrium plant biomass values under R2-limitation (equation 9)
BFbar <- (bF - uF)/(gF*uF)
B0bar <- (b0 - u0)/(g0*u0)
BFbar;B0bar

#Calculate conditions for Hopf bifurcation (Tr*SPM<Det; complex conjugate eigenvalues)

#Choose unavailable N value that creates internal equilibrium
D <- 1000 #Sb<S<Sa

#load package for matrix calculations
install.packages('matrixcalc')
library('matrixcalc')

#Internal quasi-equilibria when N-fixer is R2-limited and nonfixer is N-limited
print(BF.coex <- (I+m*D-Sa)/(AbarBF*vF-AbarB0*v0)) #Equation S11
print(B0.coex <- BFbar-BF.coex) #Equation S12
#Available N is at AbarB0

#Partial derivatives for matrix evaluated at the internal quasi-equilibrium:
print(BF.BF <- -uF^2*gF*BF.coex/bF)
print(BF.B0 <- -uF^2*gF*BF.coex/bF)
print(BF.A <- 0)
print(B0.BF <- 0)
print(B0.B0 <- 0)
print(B0.A <- B0.coex*w0*v0)
print(A.BF <- uF^2*gF*BF.coex/(wF*bF)+f-uF/wF)
print(A.B0 <- uF^2*gF*BF.coex/(wF*bF)-u0/w0)
print(A.A <- -B0.coex*v0-k)

#Set up 3x3 matrix for N-fixer, nonfixer, and available N
print(m3 <- matrix(c(BF.BF,BF.B0,BF.A,B0.BF,B0.B0,B0.A,A.BF,A.B0,A.A),3,3,byrow=T))

#Calculate eigenvalues (one is negative, other two are complex conjugates, indicating a Hopf bifurcation)
m3_eigen <- eigen(m3)
m3_eigen$values

#Also look at Routh-Hurwitz criterion for stability or 3x3 matrix (Tr<0, Det<0, and Tr*SPM<Det)
#When limit cycles arise the first two conditions for stability are satisfied (Tr<0, Det<0)
#but the third is not (Tr*SPM<Det is not satisfied)

#Calculate trace (Tr<0 for stability)
print(Tr <- matrix.trace(m3))

#Calculate determinant (Det<0 for stability)
print (Det <- det(m3))

#Pricipal minors
m3_p1 <- matrix(c(B0.B0,B0.A,A.B0,A.A),2,2,byrow=T)
m3_p1
m3_p1_d <- det(m3_p1)

m3_p2 <- matrix(c(BF.BF,BF.A,A.BF,A.A),2,2,byrow=T)
m3_p2
m3_p2_d <- det(m3_p2)

m3_p3 <- matrix(c(BF.BF,BF.B0,B0.BF,B0.B0),2,2,byrow=T)
m3_p3
m3_p3_d <- det(m3_p3)

#sum of principal minors
print(SPM <- sum(m3_p1_d,m3_p2_d,m3_p3_d))

Tr*SPM < Det #Condition for stability (not satisfied; limit cycle)

#simulate successional trajectories under static N supply
sim.S2a<-lsoda(c(24000,2000,AbarB0),ts,static.D.obl,parms)
sim.S2b<-lsoda(c(27000,2500,AbarB0),ts,static.D.obl,parms)

#Figure S2A
plot(sim.S2a[,2],sim.S2a[,3],type="l",xlab=NA,ylab=NA,las=1,lwd=1,ylim=c(0,7000))
lines(sim.S2b[,2],sim.S2b[,3],lwd=1)
points(24000,2000,pch=8,col="red")
points(27000,2500,pch=8,col="red")
mtext(expression('Nonfixer biomass (kg C ha'^-1*')'),side=2,line=3.5,cex=1.5)
mtext(expression('N-fixer biomass (kg C ha'^-1*')'),side=1,line=3,cex=1.5)

#Figure S2B
plot(sim.S2a[,2]~log10(sim.S2a[,1]),type="l",xlab=NA,ylab=NA,las=1,ylim=c(0,31000),col="red",xaxt="n",lwd=1)
lines(sim.S2a[,3]~log10(sim.S2a[,1]),col="blue",lwd=1)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
legend(x=2,y=20000,legend=c(expression('N-fixer ('*italic('B'[F])*')'),expression('Nonfixer ('*italic('B'[0])*')')),col=c("red","blue"),lty=c(1,1),bty="n")
mtext(expression('Plant biomass (kg C ha'^-1*')'),side=2,line=3.5,cex=1.5)
mtext("Time (y)",side=1,line=3,cex=1.5)


################################################################################
#Simulations in figure S3 ######################################################
################################################################################

#Coexistence stable at long-term equilibrium (Sb<S.coex<Sa)
#Figure S3 is the same simulation as figure 3E,F
#Parameters for figure S3 (changes parameters back to base values; Table S1)
I <- 4
uF <- 0.2
wF <- 50
w0 <- 50
vF <- 0.05

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#solve for critical S values that were changed for figure 3C,D and S2 (SF<Sb<Sa<SF.star<S0)
print(Sa <- (k + v0*(bF - uF)/(gF*uF))*u0/(w0*v0)) #Equation 11
print(Sb <- min((k + vF*(bF - uF)/(gF*uF))*u0/(w0*v0),k*u0/(w0*v0)+uF*((bF - uF)/(gF*uF))/wF)) #Equation 12
print(SF <- (k + vF*(bF - uF)/(gF*uF))*(uF - wF*f)/(wF*vF)) #Equation 13
print(SF.star <- (k + vF*(bF - uF)/(gF*uF))*(uF)/(wF*vF)) #Equation 13 with downregulated SNF (facultative SNF)
print(S0 <- (k + v0*(b0 - u0)/(g0*u0))*u0/(w0*v0)) #Equation 14

#equilibrium A values under N-limitation (equation 8)
AbarBF <- (uF - wF*f)/(vF*wF)
AbarB0 <- u0/(v0*w0) 
AbarBF;AbarB0

#equilibrium plant biomass values under R2-limitation (equation 9)
BFbar <- (bF - uF)/(gF*uF)
B0bar <- (b0 - u0)/(g0*u0)
BFbar;B0bar

#Additional equilibria
print(S.coex <- I+m*(BFbar*u0/w0+(uF/wF-u0/w0)*(I-Sa)/(AbarB0*(vF-v0)))/ #Equation S23
                 (m+phi-(uF/wF-u0/w0)*m/(AbarB0*(vF-v0))))
print(BF.coex <- (S.coex-Sa)/(AbarB0*(vF-v0))) #Equation S15
print(B0.coex <- BFbar-BF.coex) #Equation S16
print(AbarBF.colim <- S.coex/(k+vF*BFbar)) #Equation S7

#Figure S3A
par(pty="s",mar=c(5.1, 5.1, 4.1, 2.1))
plot(log10(sim.3E.dynamic[,1]),sim.3E.dynamic[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.3E.dynamic[,1]),sim.3E.dynamic[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
abline(h=BF.coex,lty=2,col="red")
abline(h=B0.coex,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(B)['F,coex,B'['F,colim']])),side=4,line=0.2,at=BF.coex,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)['0,coex,B'['F,colim']])),side=4,line=0.2,at=B0.coex,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S3B
par(pty="s",mar=c(5.1, 5.1, 4.1, 2.1))
plot(log10(sim.3F.dynamic[,1]),sim.3F.dynamic[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.3F.dynamic[,1]),sim.3F.dynamic[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
abline(h=BF.coex,lty=2,col="red")
abline(h=B0.coex,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(B)['F,coex,B'['F,colim']])),side=4,line=0.2,at=BF.coex,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)['0,coex,B'['F,colim']])),side=4,line=0.2,at=B0.coex,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S3C
plot(log10(sim.3E.dynamic[,1]),I+m*sim.3E.dynamic[,4],type="l",ylim=c(0,160),las=1,ylab=NA,xlab=NA,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.coex,lty=2,col="purple")
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(S[b])),side=4,line=0.2,at=Sb,col="gray50",cex=0.8,las=1)
mtext(text=expression(italic(bar(S)['coex,B'['F,colim']])),side=4,line=0.2,at=S.coex,col="purple",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S3D
plot(log10(sim.3F.dynamic[,1]),I+m*sim.3F.dynamic[,4],type="l",ylim=c(0,160),las=1,ylab=NA,xlab=NA,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.coex,lty=2,col="purple")
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(S[b])),side=4,line=0.2,at=Sb,col="gray50",cex=0.8,las=1)
mtext(text=expression(italic(bar(S)['coex,B'['F,colim']])),side=4,line=0.2,at=S.coex,col="purple",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S3E
plot(log10(sim.3E.dynamic[,1]),sim.3E.dynamic[,5],type="l",ylim=c(0,.1),ylab=NA,xlab=NA,las=1,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S3F
plot(log10(sim.3F.dynamic[,1]),sim.3F.dynamic[,5],type="l",ylim=c(0,.1),ylab=NA,xlab=NA,las=1,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

################################################################################
#Simulations in figure S5 ######################################################
################################################################################

#Figure S5 simulated before figure S4 because figure S3 and S5 both show initial coexistence (Sb<S<Sa)

#Transient coexistence transitioning into nonfixers excluding N-fixers at long-term equilibrium (figure S5;Sb<Sa<S.coex<S.B0.Nlim<S0)
#new parameters for figure S5
I <- 8
phi <- 0.0005

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#equilibrium A values under N-limitation remain unchanged (equation 8)

#equilibrium plant biomass values under R2-limitation remain unchanged (equation 9)

#Critical N supply points remain unchanged

#Equilibria
print(S.B0.Nlim <- I+m*(I-k*AbarB0)/phi) #Equation S19
print(B0.Nlim <- (w0*v0*S.B0.Nlim-k*u0)/(v0*u0)) #Equation S4
print(AbarBF.colim <- S.B0.Nlim/(k+vF*BFbar)) #Equation S7
print(S.coex <- I+m*(BFbar*u0/w0+(uF/wF-u0/w0)*(I-Sa)/(AbarB0*(vF-v0)))/ #Equation S23
        (m+phi-(uF/wF-u0/w0)*m/(AbarB0*(vF-v0))))

#simulate successional trajectories
D <- 13000 #Sb<S<Sa
sim.S5A <- lsoda(c(8000,2000,D,0),ts,dynamic.D.fac,parms)
sim.S5B <- lsoda(c(2000,8000,D,AbarB0),ts,dynamic.D.fac,parms)

#Figure S5A
plot(log10(sim.S5A[,1]),sim.S5A[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.S5A[,1]),sim.S5A[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
abline(h=B0.Nlim,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(B)['0,Nlim'])),side=4,line=0.2,at=B0.Nlim,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S5B
plot(log10(sim.S5B[,1]),sim.S5B[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.S5B[,1]),sim.S5B[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
abline(h=B0.Nlim,lty=2,col="blue")
mtext(text=expression(italic(bar(B)['0,Nlim'])),side=4,line=0.2,at=B0.Nlim,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S5C
plot(log10(sim.S5A[,1]),I+m*sim.S5A[,4],type="l",ylim=c(0,160),las=1,ylab=NA,xlab=NA,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.B0.Nlim,lty=2,col="blue")
abline(h=S.coex,lty=2,col="purple")
mtext(text=expression(italic(bar(S)[B['0,Nlim']])),side=4,line=0.2,at=S.B0.Nlim,col="blue",cex=.8,las=1)
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(S[b])),side=4,line=0.2,at=Sb,col="gray50",cex=0.8,las=1)
mtext(text=expression(italic(bar(S)['coex,B'['F,colim']])),side=4,line=0.2,at=S.coex,col="purple",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S5D
plot(log10(sim.S5B[,1]),I+m*sim.S5B[,4],type="l",ylim=c(0,160),las=1,ylab=NA,xlab=NA,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.B0.Nlim,lty=2,col="blue")
abline(h=S.coex,lty=2,col="purple")
mtext(text=expression(italic(bar(S)[B['0,Nlim']])),side=4,line=0.2,at=S.B0.Nlim,col="blue",cex=.8,las=1)
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(S[b])),side=4,line=0.2,at=Sb,col="gray50",cex=0.8,las=1)
mtext(text=expression(italic(bar(S)['coex,B'['F,colim']])),side=4,line=0.2,at=S.coex,col="purple",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S5E
plot(log10(sim.S5A[,1]),sim.S5A[,5],type="l",ylim=c(0,.1),ylab=NA,xlab=NA,las=1,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S5F
plot(log10(sim.S5B[,1]),sim.S5B[,5],type="l",ylim=c(0,0.1),xaxt="n",xlab=NA,ylab=NA,las=1,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

################################################################################
#Simulations in figure S4 ######################################################
################################################################################

#Priority effects stable at long-term equilibrium (figure S4C: Sa<S.BFbar<Sb; figure S4D: Sa<S.B0.Nlim<S0<Sb)
#Figure S4 is the same simulation as figure 3C,D
#Parameters for figure S4 (changes transient coexistence parameters [figure S5] back to figures 3C,3D,S4 parameters)
wF <- 40
vF <- 0.1

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#solve for critical S values that were changed for figure S5 (SF<Sa<S0<SF.star<Sb;Sa & S0 remain unchanged)
print(Sb <- min((k + vF*(bF - uF)/(gF*uF))*u0/(w0*v0),k*u0/(w0*v0)+uF*((bF - uF)/(gF*uF))/wF)) #Equation 12
print(SF <- (k + vF*(bF - uF)/(gF*uF))*(uF - wF*f)/(wF*vF)) #Equation 13
print(SF.star <- (k + vF*(bF - uF)/(gF*uF))*(uF)/(wF*vF)) #Equation 13 with downregulated SNF (facultative SNF)

#Equilibria
print(S.B0.Nlim <- I+m*(I-k*AbarB0)/phi) #Equation S19
print(B0.Nlim <- (w0*v0*S.B0.Nlim-k*u0)/(v0*u0)) #Equation S4
print(S.BFbar <- I+m*uF*BFbar/(wF*(phi+m))) #Equation S18
print(AbarBF.colim.A <- S.BFbar/(k+vF*BFbar)) #Equation S7; figure S4E
print(AbarBF.colim.B <- S.B0.Nlim/(k+vF*BFbar)) #Equation S7; figure S4F

#Figure S4A
plot(log10(sim.3C.dynamic[,1]),sim.3C.dynamic[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.3C.dynamic[,1]),sim.3C.dynamic[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S4B
plot(log10(sim.3D.dynamic[,1]),sim.3D.dynamic[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.3D.dynamic[,1]),sim.3D.dynamic[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
abline(h=B0.Nlim,lty=2,col="blue")
mtext(text=expression(italic(bar(B)['0,Nlim'])),side=4,line=0.2,at=B0.Nlim,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S4C
plot(log10(sim.3C.dynamic[,1]),I+m*sim.3C.dynamic[,4],type="l",ylim=c(0,200),las=1,ylab=NA,xlab=NA,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.BFbar,col="red",lty=2)
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic('S'['b'])),side=4,line=0.2,at=Sb,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(bar(S)[bar(B)[F]])),side=4,line=0.2,at=S.BFbar,col="red",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S4D
plot(log10(sim.3D.dynamic[,1]),I+m*sim.3D.dynamic[,4],type="l",ylim=c(0,200),las=1,xaxt="n",ylab=NA,xlab=NA,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.B0.Nlim,lty=2,col="blue")
mtext(text=expression(italic(bar(S)[B['0,Nlim']])),side=4,line=0.2,at=S.B0.Nlim,col="blue",cex=.8,las=1)
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic('S'['b'])),side=4,line=0.2,at=Sb,col="gray50",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S4E
plot(log10(sim.3C.dynamic[,1]),sim.3C.dynamic[,5],type="l",ylim=c(0,.1),ylab=NA,xlab=NA,las=1,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim.A,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim.A,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S4F
plot(log10(sim.3D.dynamic[,1]),sim.3D.dynamic[,5],type="l",ylim=c(0,0.1),xaxt="n",xlab=NA,ylab=NA,las=1,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim.B,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim.B,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

################################################################################
#Simulations in figure S10 ######################################################
################################################################################

#Figure S10 simulated before figures S6-S8 because it is based on figure S4

#Equation 7
sim.Ff <- rep(NA,length(sim.3C.dynamic[,1]))
for(i in 1:length(sim.3C.dynamic[,1])){
  sim.Ff[i] <- max(0,min(bF/(wF*(1 + gF*(sim.3C.dynamic[i,2] + sim.3C.dynamic[i,3]))) - vF*sim.3C.dynamic[i,5],f))
}

#Figure S10A
par(pty="s",mar=c(5.1, 8.1, 4.1, 2.1))
plot(sim.Ff~log10(sim.3C.dynamic[,1]),type="l",xaxt="n",las=1,xlab=NA,ylab=NA,lwd=2,ylim=c(0,0.01))
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
mtext(expression('N-fixation'),side=2,line=7,cex=1.5)
mtext(expression('per fixer biomass'),side=2,line=5,cex=1.5)
mtext(expression('('*italic(F)*'; kg N kg C'^-1*' y'^-1*')'),side=2,line=3.5,cex=1.2)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S10B
plot(sim.Ff*sim.3C.dynamic[,2]~log10(sim.3C.dynamic[,1]),type="l",xaxt="n",las=1,xlab=NA,ylab=NA,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
mtext(expression('N-fixation'),side=2,line=5,cex=1.5)
mtext(expression('('*italic(B[F])*''%*%''*italic(F)*'; kg N ha'^-1*' y'^-1*')'),side=2,line=3,cex=1.2)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S10C
plot((I+m*sim.3C.dynamic[,4])/(k+vF*BFbar)~log10(sim.3C.dynamic[,1]),type="l",xaxt="n",las=1,xlab=NA,ylab=NA,lwd=2,ylim=c(0,0.1))
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarBF.colim.A,lty=2,col="red")
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim.A,col="red",cex=.8,las=1)
mtext(expression('Quasi-equilibrium'),side=2,line=7,cex=1.5)
mtext(expression('of available N'),side=2,line=5.5,cex=1.5)
mtext(expression('('*italic(hat(A)[B['F,colim']])*'; kg N ha'^-1*')'),side=2,line=3,cex=1.2)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

################################################################################
#Simulations in figure S6 ######################################################
################################################################################

#Transient priority effects transitioning into N-fixers excluding nonfixers at long-term equilibrium (figure S6C: Sa<S.BFbar<Sb; figure S6D: S.B0.Nlim<Sa<S.BFbar<Sb)
#new parameter for figure S6
I <- 2

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#equilibrium A values under N-limitation remain unchanged (equation 8)

#equilibrium plant biomass values under R2-limitation remain unchanged (equation 9)

#Critical N supply points remain unchanged

#Equilibria
print(S.B0.Nlim <- I+m*(I-k*AbarB0)/phi) #Equation S19
print(B0.Nlim <- (w0*v0*S.B0.Nlim-k*u0)/(v0*u0)) #Equation S4
print(S.BFbar <- I+m*uF*BFbar/(wF*(phi+m))) #Equation S18
print(AbarBF.colim <- S.BFbar/(k+vF*BFbar)) #Equation S7

#simulate successional trajectories
D <- 16000 #Sa<S<Sb
sim.S6A <- lsoda(c(8000,2000,D,0),ts,dynamic.D.fac,parms)
sim.S6B <- lsoda(c(2000,8000,D,AbarB0),ts,dynamic.D.fac,parms)

#Figure S6A
par(pty="s",mar=c(5.1, 5.1, 4.1, 2.1))
plot(log10(sim.S6A[,1]),sim.S6A[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.S6A[,1]),sim.S6A[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S6B
plot(log10(sim.S6B[,1]),sim.S6B[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.S6B[,1]),sim.S6B[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S6C
plot(log10(sim.S6A[,1]),I+m*sim.S6A[,4],type="l",ylim=c(0,200),las=1,ylab=NA,xlab=NA,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.BFbar,col="red",lty=2)
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic('S'['b'])),side=4,line=0.2,at=Sb,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(bar(S)[bar(B)[F]])),side=4,line=0.2,at=S.BFbar,col="red",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S6D
plot(log10(sim.S6B[,1]),I+m*sim.S6B[,4],type="l",ylim=c(0,200),las=1,xaxt="n",ylab=NA,xlab=NA,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.B0.Nlim,lty=2,col="blue")
abline(h=S.BFbar,col="red",lty=2)
mtext(text=expression(italic(bar(S)[B['0,Nlim']])),side=4,line=0.2,at=S.B0.Nlim,col="blue",cex=.8,las=1)
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic('S'['b'])),side=4,line=0.2,at=Sb,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(bar(S)[bar(B)[F]])),side=4,line=0.2,at=S.BFbar,col="red",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S6E
plot(log10(sim.S6A[,1]),sim.S6A[,5],type="l",ylim=c(0,.1),ylab=NA,xlab=NA,las=1,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S6F
plot(log10(sim.S6B[,1]),sim.S6B[,5],type="l",ylim=c(0,0.1),xaxt="n",xlab=NA,ylab=NA,las=1,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

################################################################################
#Simulations in figure S7 ######################################################
################################################################################

#Transient priority effects transitioning into nonfixers excluding N-fixers at long-term equilibrium (figure S7C:Sa<S.B0bar<Sb<S.BFbar<S0<S.B0.Nlim; figure S7D: Sa<S.B0bar<Sb<S0<S.B0.Nlim)
#new parameter for figure S7
I <- 14

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#equilibrium A values under N-limitation remain unchanged (equation 8)

#equilibrium plant biomass values under R2-limitation remain unchanged (equation 9)

#Critical N supply points remain unchanged

#Equilibria
print(S.BFbar <- I+m*uF*BFbar/(wF*(phi+m))) #Equation S18
print(S.B0bar <- I+m*u0*B0bar/(w0*(phi+m))) #Equation S20
print(AbarBF.colim <- S.B0bar/(k+vF*BFbar)) #Equation S7
print(AbarB0.R2lim <- (S.B0bar-u0*B0bar/w0)/k) #Equation S5

#simulate successional trajectories
D <- 16000 #Sa<S<Sb
sim.S7A <- lsoda(c(8000,2000,D,0),ts,dynamic.D.fac,parms)
sim.S7B <- lsoda(c(2000,8000,D,AbarB0),ts,dynamic.D.fac,parms)

#Figure S7A
plot(log10(sim.S7A[,1]),sim.S7A[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.S7A[,1]),sim.S7A[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S7B
plot(log10(sim.S7B[,1]),sim.S7B[,2],type="l",col="red",lwd=2,ylab=NA,xlab=NA,ylim=c(0,45000),las=1,xaxt="n")
lines(log10(sim.S7B[,1]),sim.S7B[,3],col="blue",lty=2,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=BFbar,lty=2,col="red")
abline(h=B0bar,lty=2,col="blue")
mtext(text=expression(italic(bar(B)[F])),side=4,line=0.2,at=BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(B)[0])),side=4,line=0.2,at=B0bar,col="blue",cex=.8,las=1)
mtext(expression('Plant Biomass'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg C ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S7C
plot(log10(sim.S7A[,1]),I+m*sim.S7A[,4],type="l",ylim=c(0,200),las=1,ylab=NA,xlab=NA,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.BFbar,col="red",lty=2)
abline(h=S.B0bar,col="blue",lty=2)
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic('S'['b'])),side=4,line=0.2,at=Sb,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(bar(S)[bar(B)[F]])),side=4,line=0.2,at=S.BFbar,col="red",cex=.8,las=1)
mtext(text=expression(italic(bar(S)[bar(B)[0]])),side=4,line=0.2,at=S.B0bar,col="blue",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S7D
plot(log10(sim.S7B[,1]),I+m*sim.S7B[,4],type="l",ylim=c(0,200),las=1,xaxt="n",ylab=NA,xlab=NA,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=Sa,col="gray50")
abline(h=Sb,col="gray50")
abline(h=S.B0bar,lty=2,col="blue")
mtext(text=expression(italic('S'['a'])),side=4,line=0.2,at=Sa,col="gray50",cex=.8,las=1)
mtext(text=expression(italic('S'['b'])),side=4,line=0.2,at=Sb,col="gray50",cex=.8,las=1)
mtext(text=expression(italic(bar(S)[bar(B)[0]])),side=4,line=0.2,at=S.B0bar,col="blue",cex=.8,las=1)
mtext(expression('N supply'),side=2,line=5.5,cex=1.5)
mtext(expression('(kg N ha'^-1*' y'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S7E
plot(log10(sim.S7A[,1]),sim.S7A[,5],type="l",ylim=c(0,.1),ylab=NA,xlab=NA,las=1,xaxt="n",lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

#Figure S7F
plot(log10(sim.S7B[,1]),sim.S7B[,5],type="l",ylim=c(0,0.1),xaxt="n",xlab=NA,ylab=NA,las=1,lwd=2)
axis(1,at=c(0,1,2,3,4),labels=c(1,10,100,1000,10000))
abline(h=AbarB0,col="blue",lty=2)
abline(h=AbarBF.colim,lty=2,col="red")
mtext(text=expression(italic(bar(A)['0'])),side=4,line=0.2,at=AbarB0,col="blue",cex=.8,las=1)
mtext(text=expression(italic(bar(A)[B['F,colim']])),side=4,line=0.2,at=AbarBF.colim,col="red",cex=.8,las=1)
mtext(expression('Available N'),side=2,cex=1.5,line=5.5)
mtext(expression('(kg N ha'^-1*')'),side=2,cex=1.5,line=3.5)
mtext(expression('Time (y)'),side=1,line=3,cex=1.5)

################################################################################
#Simulations in figure S9 ######################################################
################################################################################

#Set parameter I back to what it was in figures 3C,3D,S4,S10
I <- 8

parms <- list(I,m,k,uF,u0,wF,w0,vF,v0,f,phi,bF,b0,gF,g0)

#Sequences for starting N-fixer and nonfixer values
BF.seq <- seq(0,50000,5000)
B0.seq <- seq(0,50000,5000)

#Set unavailable N pool size for figure S9A (low N supply where priority effects occur)
D <- 15200 #Sa<S<Sb

#Make N-fixer and nonfixer arrays to fill with simulation across grid of starting conditions
BF.lowS <- array(NA,dim=c(10,10,length(ts)))
B0.lowS <- array(NA,dim=c(10,10,length(ts)))

#Simulate successional trajectories at static N supply and fill arrays with N-fixer and nonfixer values
for (i in 1:dim(BF.lowS)[1]){
  for (j in 1:dim(BF.lowS)[2]){
  lowS.sim <- lsoda(c(BF.seq[i],B0.seq[j],AbarB0),ts,static.D.fac,parms)
  BF.lowS[i,j,] <- lowS.sim[,2]
  B0.lowS[i,j,] <- lowS.sim[,3]
  }
}

#Matrix of colors that makes successional trajectory red if N-fixer wins and blue if nonfixer wins
col.lowS <- matrix(NA,nrow=10,ncol=10)

for (i in 1:dim(col.lowS)[1]){
  for (j in 1:dim(col.lowS)[2]){
    if (BF.lowS[i,j,length(ts)]>B0.lowS[i,j,length(ts)]){
      col.lowS[i,j] <- "red"
    } else{
      col.lowS[i,j] <- "blue"
    }
  }
}

#Set unavailable N pool size for figure S9B
D <- 17200 #Sa<S<Sb

#Make N-fixer and nonfixer arrays to fill with simulation across grid of starting conditions
BF.highS <- array(NA,dim=c(10,10,length(ts)))
B0.highS <- array(NA,dim=c(10,10,length(ts)))

#Make N-fixer and nonfixer arrays to fill with simulation across grid of starting conditions
for (i in 1:dim(BF.highS)[1]){
  for (j in 1:dim(BF.highS)[2]){
    highS.sim <- lsoda(c(BF.seq[i],B0.seq[j],AbarB0),ts,static.D.fac,parms)
    BF.highS[i,j,] <- highS.sim[,2]
    B0.highS[i,j,] <- highS.sim[,3]
  }
}

#Matrix of colors that makes successional trajectory red if N-fixer wins and blue if nonfixer wins
col.highS <- matrix(NA,nrow=10,ncol=10)

for (i in 1:dim(col.highS)[1]){
  for (j in 1:dim(col.highS)[2]){
    if (BF.highS[i,j,length(ts)]>B0.highS[i,j,length(ts)]){
      col.highS[i,j] <- "red"
    } else{
      col.highS[i,j] <- "blue"
    }
  }
}

#Figure S9A
plot(1, type="n", xlab=NA, ylab=NA, xlim=c(0,45000), ylim=c(0,45000),las=1)
for (i in 1:dim(BF.lowS)[1]){
  for (j in 1:dim(BF.lowS)[2]){
    points(BF.lowS[i,j,]~B0.lowS[i,j,],pch=20,cex=0.2,col=col.lowS[i,j])
  }
}
points(0,tail(BF.lowS[10,1,],1),pch=19,cex=1.2,col="red")
points(tail(B0.lowS[1,10,],1),0,pch=19,cex=1.2,col="blue")
points(0,0,pch=19,cex=1.2,col="white")
points(0,tail(BF.lowS[10,1,],1),pch=1,cex=1.2)
points(tail(B0.lowS[1,10,],1),0,pch=1,cex=1.2)
points(0,0,pch=1,cex=1.2)
mtext(expression('N-fixer biomass (kg C ha'^-1*')'),side=2,line=3.5,cex=1.5)
mtext(expression('Nonfixer biomass (kg C ha'^-1*')'),side=1,line=3,cex=1.5)

#Figure S9B
plot(1, type="n", xlab=NA, ylab=NA, xlim=c(0,45000), ylim=c(0,45000),las=1)
for (i in 1:dim(BF.highS)[1]){
  for (j in 1:dim(BF.highS)[2]){
    points(BF.highS[i,j,]~B0.highS[i,j,],pch=20,cex=0.2,col=col.highS[i,j])
  }
}
points(0,tail(BF.highS[10,1,],1),pch=19,cex=1.2,col="red")
points(tail(B0.highS[1,10,],1),0,pch=19,cex=1.2,col="blue")
points(0,0,pch=19,cex=1.2,col="white")
points(0,tail(BF.highS[10,1,],1),pch=1,cex=1.2)
points(tail(B0.highS[1,10,],1),0,pch=1,cex=1.2)
points(0,0,pch=1,cex=1.2)
mtext(expression('N-fixer biomass (kg C ha'^-1*')'),side=2,line=3.5,cex=1.5)
mtext(expression('Nonfixer biomass (kg C ha'^-1*')'),side=1,line=3,cex=1.5)
