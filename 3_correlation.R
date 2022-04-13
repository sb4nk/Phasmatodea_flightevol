########################################################################################

# R code used in Bank & Bradler (2022)
# https://github.com/sb4nk/Phasmatodea_flightevol

# Testing correlation between evolution of wings and ocelli.

########################################################################################

library(phytools)
library(geiger)


setwd("C:/Users/dusty/Dropbox/PhD/MS_Bank-Bradler_Ocelli-Evolution/GitHub/3_correlation/")
#setwd("~/ownCloud/Phylogenetics/Ocelli/11_ASR")

# Load tree file
tree <- read.nexus("treefile.BI_constraintsB2.nex")

### Load morphological data (wings/ocelli) for males and females
data.m <- read.csv("correlation.test.wingstates.csv", row.names=1)
x<-setNames(data.m[[3]],rownames(data.m)) # datawings
y<-setNames(data.m[[2]],rownames(data.m)) # ocelli




# NO DEPENDENT VARIABLE

fit.no<-fitPagel(tree,x,y, model="ARD",lwd.by.rate=TRUE)
plot(fit.no,lwd.by.rate=TRUE)
# only wings
w<-setNames(data.wings[[1]],rownames(data.m)) # wings
fit.wings<-fitDiscrete(tree,w, model="ARD",lwd.by.rate=TRUE)
plot(fit.wings,lwd.by.rate=TRUE)#,lwd.by.rate=TRUE)



######################## REAL STUFF


# permits (respectively) the substitution rate of x to depend on the state of y and vice versa
fit.xy<-fitPagel(tree,x,y, dep.var = "xy", model="ARD",
                 lwd.by.rate=TRUE,
                 main=c("a) Independent model","b) Full dependent model"))
fit.xy
plot(fit.xy,lwd.by.rate=TRUE)

#permits the substitution rate of x (only) to depend on the state of y
fit.x<-fitPagel(tree,x,y,dep.var="x",model="ARD",
                lwd.by.rate=TRUE,
                main=c("a) Independent model","b) Dependent model (winged | wingless)"))
fit.x
plot(fit.x,lwd.by.rate=TRUE)

#permits the substitution rate of y (only) to depend on the state of x
fit.y<-fitPagel(tree,x,y,dep.var="y",model="ARD",
                lwd.by.rate=TRUE,
                main=c("a) Independent model","b) Dependent model (ocelli | no ocelli)"))
fit.y
plot(fit.y,lwd.by.rate=TRUE)


aic<-setNames(c(fit.xy$independent.AIC,
                fit.x$dependent.AIC,
                fit.y$dependent.AIC,
                fit.xy$dependent.AIC),
              c("independent","dependent x",
                "dependent y","dependent x&y"))
aic

aic.w<-function(aic,signif=4){
  d.aic<-aic-min(aic)
  round(exp(-1/2*d.aic)/sum(exp(-1/2*d.aic)),signif)
}
aic.w(aic)


##############################################################################

# changing root state

# add root branch
rootedtree<-bind.tip(tree,tip.label="ROOT",edge.length=0,
                  where=Ntip(tree)+1)
plotTree(rootedtree,ftype="off",lwd=1)
nn<-which(rootedtree$tip.label=="ROOT")
tiplabels(rootedtree$tip.label[nn],nn,adj=c(-0.1,0.5),
          frame="none",cex=0.8,font=3)
plotTree(rootedtree,lwd=0.5)


### Root = no wings | no ocelli ###############################################

# add root state no wings - no ocelli
y.ocelli<-as.factor(c(setNames("absent","ROOT"),
                      setNames(as.character(y),names(y))))
x.wings<-as.factor(c(setNames("0","ROOT"),
                     setNames(as.character(x),names(x))))

# permits (respectively) the substitution rate of x to depend on the state of y and vice versa
fit.xy1<-fitPagel(rootedtree,x.wings,y.ocelli, dep.var = "xy", model="ARD",
                 lwd.by.rate=TRUE,
                 main=c("a) Independent model","b) Full dependent model"))
fit.xy1
plot(fit.xy1,lwd.by.rate=TRUE)

#permits the substitution rate of x (only) to depend on the state of y
fit.x1<-fitPagel(rootedtree,x.wings,y.ocelli,dep.var="x",model="ARD",
                lwd.by.rate=TRUE,
                main=c("a) Independent model","b) Dependent model (winged | wingless)"))
fit.x1
plot(fit.x1,lwd.by.rate=TRUE)

#permits the substitution rate of y (only) to depend on the state of x
fit.y1<-fitPagel(rootedtree,x.wings,y.ocelli,dep.var="y",model="ARD",
                lwd.by.rate=TRUE,
                main=c("a) Independent model","b) Dependent model (ocelli | no ocelli)"))
fit.y1
plot(fit.y1,lwd.by.rate=TRUE)


aic<-setNames(c(fit.xy1$independent.AIC,
                fit.x1$dependent.AIC,
                fit.y1$dependent.AIC,
                fit.xy1$dependent.AIC),
              c("independent","dependent x",
                "dependent y","dependent x&y"))
aic

aic.w<-function(aic,signif=4){
  d.aic<-aic-min(aic)
  round(exp(-1/2*d.aic)/sum(exp(-1/2*d.aic)),signif)
}
aic.w(aic)



### Root = wings | no ocelli ###############################################

# add root state wings - no ocelli
y.ocelli<-as.factor(c(setNames("absent","ROOT"),
                      setNames(as.character(y),names(y))))
x.wings<-as.factor(c(setNames("1","ROOT"),
                     setNames(as.character(x),names(x))))

# permits (respectively) the substitution rate of x to depend on the state of y and vice versa
fit.xy2<-fitPagel(rootedtree,x.wings,y.ocelli, dep.var = "xy", model="ARD",
                  lwd.by.rate=TRUE,
                  main=c("a) Independent model","b) Full dependent model"))
fit.xy2
plot(fit.xy2,lwd.by.rate=TRUE)

#permits the substitution rate of x (only) to depend on the state of y
fit.x2<-fitPagel(rootedtree,x.wings,y.ocelli,dep.var="x",model="ARD",
                 lwd.by.rate=TRUE,
                 main=c("a) Independent model","b) Dependent model (winged | wingless)"))
fit.x2
plot(fit.x2,lwd.by.rate=TRUE)

#permits the substitution rate of y (only) to depend on the state of x
fit.y2<-fitPagel(rootedtree,x.wings,y.ocelli,dep.var="y",model="ARD",
                 lwd.by.rate=TRUE,
                 main=c("a) Independent model","b) Dependent model (ocelli | no ocelli)"))
fit.y2
plot(fit.y2,lwd.by.rate=TRUE)


aic<-setNames(c(fit.xy2$independent.AIC,
                 fit.x2$dependent.AIC,
                 fit.y2$dependent.AIC,
                 fit.xy2$dependent.AIC),
               c("independent","dependent x",
                 "dependent y","dependent x&y"))
aic

aic.w<-function(aic,signif=4){
  d.aic<-aic-min(aic)
  round(exp(-1/2*d.aic)/sum(exp(-1/2*d.aic)),signif)
}
aic.w(aic)



### Root = wings | ocelli ###############################################

# add root state wings - ocelli
y.ocelli<-as.factor(c(setNames("present","ROOT"),
                      setNames(as.character(y),names(y))))
x.wings<-as.factor(c(setNames("1","ROOT"),
                     setNames(as.character(x),names(x))))

# permits (respectively) the substitution rate of x to depend on the state of y and vice versa
fit.xy3<-fitPagel(rootedtree,x.wings,y.ocelli, dep.var = "xy", model="ARD",
                  lwd.by.rate=TRUE,
                  main=c("a) Independent model","b) Full dependent model"))
fit.xy3
plot(fit.xy2,lwd.by.rate=TRUE)

#permits the substitution rate of x (only) to depend on the state of y
fit.x3<-fitPagel(rootedtree,x.wings,y.ocelli,dep.var="x",model="ARD",
                 lwd.by.rate=TRUE,
                 main=c("a) Independent model","b) Dependent model (winged | wingless)"))
fit.x3
plot(fit.x3,lwd.by.rate=TRUE)

#permits the substitution rate of y (only) to depend on the state of x
fit.y3<-fitPagel(rootedtree,x.wings,y.ocelli,dep.var="y",model="ARD",
                 lwd.by.rate=TRUE,
                 main=c("a) Independent model","b) Dependent model (ocelli | no ocelli)"))
fit.y3
plot(fit.y3,lwd.by.rate=TRUE)


aic<-setNames(c(fit.xy3$independent.AIC,
                fit.x3$dependent.AIC,
                fit.y3$dependent.AIC,
                fit.x3$dependent.AIC),
              c("independent","dependent x",
                "dependent y","dependent x&y"))
aic

aic.w<-function(aic,signif=4){
  d.aic<-aic-min(aic)
  round(exp(-1/2*d.aic)/sum(exp(-1/2*d.aic)),signif)
}
aic.w(aic)





########################################################################################
########################################################################################


# Testing correlation of ocelli and either of the three wing states (wingless, 
# partially-winged, fully-winged)


data.m <- read.csv("correlation.test.wingstates.csv", row.names=1)
x1<-setNames(data.m[[1]],rownames(data.m))  ## wingless
x2<-setNames(data.m[[2]],rownames(data.m))  ## partially-winged
x3<-setNames(data.m[[3]],rownames(data.m))  ## fully-winged
yy<-setNames(data.m[[4]],rownames(data.m))  ## dependent variable: ocelli



#permits the substitution rate of y (only) to depend on the state of x
fit.yy1<-fitPagel(tree,x1,yy,dep.var="y",model="ARD",
                lwd.by.rate=TRUE, method = "fitMk",
                main=c("a) Independent model","b) Dependent model (ocelli | no ocelli)"))
fit.yy1
plot(fit.yy1,lwd.by.rate=TRUE)

fit.yy2<-fitPagel(tree,x2,yy,dep.var="y",model="ARD",
                 lwd.by.rate=TRUE, method = "fitMk",
                 main=c("a) Independent model","b) Dependent model (ocelli | no ocelli)"))
fit.yy2
plot(fit.yy2,lwd.by.rate=TRUE)

fit.yy3<-fitPagel(tree,x3,yy,dep.var="y",model="ARD",
                 lwd.by.rate=TRUE, method = "fitMk",
                 main=c("a) Independent model","b) Dependent model (ocelli | no ocelli)"))
fit.yy3
plot(fit.yy3,lwd.by.rate=TRUE)


# bonferroni corrected p-values for multiple comparisons
pvalue <- 0.319379
corrected.p<-p.adjust(pvalue, method = "bonferroni", n = 3)
corrected.p 


