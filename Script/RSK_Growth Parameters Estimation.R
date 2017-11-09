setwd("~/Desktop/RSK_Script")
getwd()
rm()
ls()

library(FSA)
library(nlstools)
library(Rcmdr)

age <- read.csv("RSK.csv", as.is=TRUE)

age$Sex <- as.factor(age$Sex)
male <- subset(age, Sex=="M")
female <- subset(age, Sex=="F")

##############
### VBGM-2 ###
##############

# L0=90 mm

vbL0 <- function(t,Linf,K=NULL,L0=90) {
  if (length(Linf)==2) { 
    K <- Linf[[2]]
    Linf <- Linf[[1]]
  }
  Linf-(Linf-L0)*exp(-K*t)
}
vbL0(0,Linf=620,K=0.2)
vbL0(0,Linf=c(620,0.2))

### Female ###

svL0_f <- list(Linf=max(female$TL,na.rm=TRUE), K=0.2)
fitL0_f <- nls(TL~vbL0(Age,Linf,K),data=female,start=svL0_f)
coef(fitL0_f)

boot_f <- nlsBoot(fitL0_f)
x <- seq(0,10,length.out=199)
LCI <- UCI <- numeric(length(x))
for(i in 1:length(x)) {
  tmp <- apply(boot_f$coefboot, MARGIN=1, FUN=vbL0, t=x[i])
  LCI[i] <- quantile(tmp,0.025)
  UCI[i] <- quantile(tmp,0.975)
}

#par(mfrow=c(1,2),mar=c(2,2,0,0),oma=c(0,0,0,0))
#par(mfrow=c(2,2),mar=c(2,2,0,0),oma=c(0,0,0,0))

par(mfrow=c(2,2),mar=c(4,4,0,0),oma=c(0,0,0,0))
mtext("Age (years)", side=1, outer = FALSE)

pL0_f <- vbL0(x,Linf= coef(fitL0_f))
xlmts <- range(c(x,female$Age))
ylmts <- range(c(pL0_f,LCI,UCI,female$TL))
plot(TL~Age, data=female, xlab="", bty="l", ylab="",
     xlim=xlmts, ylim=c(0,600),pch=19,col=rgb(0,0,0,1/3))
axis(1, at = seq(0, 10, by = 1), las=1)
lines(pL0_f~x,lwd=2)
lines(UCI~x,lwd=2,lty="dashed")
lines(LCI~x,lwd=2,lty="dashed")
legend("bottomright", c("Female","95% CI"), bty = "n", border = "white", lty=c(1,2), pch=c(19))

### Male ###

svL0_m <- list(Linf=max(male$TL,na.rm=TRUE), K=0.2)
fitL0_m <- nls(TL~vbL0(Age,Linf,K),data=male,start=svL0_m)
coef(fitL0_m)

boot_m <- nlsBoot(fitL0_m)
x <- seq(0,6,length.out=199)
LCI <- UCI <- numeric(length(x))
for(i in 1:length(x)) {
  tmp <- apply(boot_m$coefboot, MARGIN=1, FUN=vbL0, t=x[i])
  LCI[i] <- quantile(tmp,0.025)
  UCI[i] <- quantile(tmp,0.975)
}

pL0_m <- vbL0(x,Linf=coef(fitL0_m))
xlmts <- range(c(x,male$Age))
ylmts <- range(c(pL0_m,LCI,UCI,male$TL))
plot(TL~Age, data=male, xlab="", bty="l", ylab="",
     xlim=xlmts, ylim=c(0,600),pch=19,col=rgb(0,0,0,1/3))
lines(pL0_m~x,lwd=2)
lines(UCI~x,lwd=2,lty="dashed")
lines(LCI~x,lwd=2,lty="dashed")
legend("bottomright", c("Male","95% CI"), bty = "n", border = "white", lty=c(1,2), pch=c(19))
legend("topright", c("(a)"),bty = "n", border = "white")

#par(mfrow=c(1,1))

##############
### VBGM-3 ###
##############
vbTyp <-vbFuns()
vbTyp(3, Linf=620,K=0.2,t0=0)

### Female ###

svTyp_f <- list(Linf=max(female$TL,na.rm=TRUE), K=0.2, t0=0)
fitTyp_f <- nls(TL~vbTyp(Age,Linf,K,t0),data=female,start=svTyp_f)

boot_f <- nlsBoot(fitTyp_f)
LCI <- UCI <- numeric(length(Age))
x <- seq(0,10,length.out=199)
for(i in 1:length(x)) {
  tmp <- apply(boot_f$coefboot, MARGIN=1, FUN=vbTyp, t=x[i])
  LCI[i] <- quantile(tmp,0.025)
  UCI[i] <- quantile(tmp,0.975)
}

#par(mfrow=c(1,2))
#par(mar=c(2,2,0,0))
#par(oma=c(0,0,0,0))

pTyp_f <- vbTyp(x,Linf=coef(fitTyp_f))
xlmts <- range(c(x,female$Age))
ylmts <- range(c(pTyp_f,LCI,UCI,female$TL))
plot(TL~Age, data=female, xlab="", bty="l", ylab="",
     xlim=xlmts, ylim=c(0,600),pch=19,col=rgb(0,0,0,1/3))
axis(1, at = seq(0, 10, by = 1), las=1)
lines(pTyp_f~x,lwd=2)
lines(UCI~x,lwd=2,lty="dashed")
lines(LCI~x,lwd=2,lty="dashed")
legend("bottomright", c("Female","95% CI"), bty = "n", border = "white", lty=c(1,2), pch=c(19))

### Male ###

svTyp_m <- list(Linf=max(male$TL,na.rm=TRUE), K=0.2, t0=0)
fitTyp_m <- nls(TL~vbTyp(Age,Linf,K,t0),data=male,start=svTyp_m)

boot_m <- nlsBoot(fitTyp_m)
LCI <- UCI <- numeric(length(Age))
x <- seq(0,6,length.out=199)
for(i in 1:length(x)) {
  tmp <- apply(boot_m$coefboot, MARGIN=1, FUN=vbTyp, t=x[i])
  LCI[i] <- quantile(tmp,0.025)
  UCI[i] <- quantile(tmp,0.975)
}

pTyp_m <- vbTyp(x,Linf=coef(fitTyp_m))
xlmts <- range(c(x,male$Age))
ylmts <- range(c(pTyp_m,LCI,UCI,male$TL))
plot(TL~Age, data=male, xlab="", bty="l", ylab="",
     xlim=xlmts, ylim=c(0,600),pch=19,col=rgb(0,0,0,1/3))
lines(pTyp_m~x,lwd=2)
lines(UCI~x,lwd=2,lty="dashed")
lines(LCI~x,lwd=2,lty="dashed")
legend("bottomright", c("Male","95% CI"), bty = "n", border = "white", lty=c(1,2), pch=c(19))
legend("topright", c("(b)"),bty = "n", border = "white")

#mtext('Age (years)', side = 1, outer = FALSE, line = 0)
#mtext('Total Length (mm)', side = 2, outer = FALSE, line = 2)

#mtext("Age (years)", side=1, outer = FALSE)

#par(mfrow=c(1,1))
