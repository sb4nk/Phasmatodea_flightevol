########################################################################################

# R code used in Bank & Bradler (2022)
# https://github.com/sb4nk/Phasmatodea_flightevol

# Ancestral states reconstruction including model testing and stochastic character 
# mapping for multi-state and binary datasets, as well as the evaluation of the number
# of transitions events.

########################################################################################

library("phytools")
library("geiger")


# Load tree file
tree <- read.nexus("treefile.BI_constraintsB2.nex")

### Load morphological data (wings/ocelli) for males and females
wing.table.m <- read.csv("morphdata_males.txt", row.names = 1)
wing.table.f <- read.csv("morphdata_females.txt", row.names = 1)

# wings
wing.m<-setNames(wing.table.m[[1]],rownames(wing.table.m))
wing.f<-setNames(wing.table.f[[1]],rownames(wing.table.f))

# ocelli
oc.m<-setNames(wing.table.m[[2]],rownames(wing.table.m))
oc.f<-setNames(wing.table.f[[2]],rownames(wing.table.f))


# on Linux systen: factorise
wing.m <- as.factor(wing.m)
wing.f <- as.factor(wing.f)
oc.m <- as.factor(oc.m)
oc.f <- as.factor(oc.f)


# include uncertainty for unknown wing states
# MALE WINGS
    wing.m.matrix.unc <- to.matrix(wing.m[tree$tip.label],levels(wing.m))
    wing.m.matrix.unc
    
    wing.m.matrix.unc["PHA4",c(1,3)]<-c(0.5,0.5)
    wing.m.matrix.unc["PHA4",2]<-0
    wing.m.matrix.unc["SB0362",]<-rep(1/3,3)
    wing.m.matrix.unc["SB0540",]<-rep(1/3,3)
    wing.m.matrix.unc["SB0481",c(1,3)]<-c(0.5,0.5)
    wing.m.matrix.unc["SB0481",2]<-0
    wing.m.matrix.unc["SB0521",]<-rep(1/3,3)
    wing.m.matrix.unc["SB0544",]<-rep(1/3,3)
    wing.m.matrix.unc["SB0571",c(1,3)]<-c(0.5,0.5)
    wing.m.matrix.unc["SB0571",2]<-0
    wing.m.matrix.unc["GB103",]<-rep(1/3,3)
    wing.m.matrix.unc["GB153",]<-rep(1/3,3)

# FEMALE WINGS
    wing.f.matrix.unc <- to.matrix(wing.f[tree$tip.label],levels(wing.f))
    wing.f.matrix.unc
    
    wing.f.matrix.unc["SB0280",]<-rep(1/3,3)
    wing.f.matrix.unc["SB0470",]<-rep(1/3,3)
    wing.f.matrix.unc["SB0476",]<-rep(1/3,3)
    wing.f.matrix.unc["SB0480",]<-rep(1/3,3)
    wing.f.matrix.unc["SB0529",]<-rep(1/3,3)
    wing.f.matrix.unc["SB0530",]<-rep(1/3,3)
    wing.f.matrix.unc["SB0536",]<-rep(1/3,3)
    wing.f.matrix.unc["SB0559",]<-rep(1/3,3)
    wing.f.matrix.unc["GB111",]<-rep(1/3,3)

# MALE OCELLI
    oc.m.matrix.unc <- to.matrix(oc.m[tree$tip.label],levels(oc.m))
    oc.m.matrix.unc
    
    oc.m.matrix.unc["PHA4",1:2]<-c(0.5,0.5)
    oc.m.matrix.unc["SB0362",1:2]<-c(0.5,0.5)
    oc.m.matrix.unc["SB0481",1:2]<-c(0.5,0.5)
    oc.m.matrix.unc["SB0521",1:2]<-c(0.5,0.5)
    oc.m.matrix.unc["SB0544",1:2]<-c(0.5,0.5)
    oc.m.matrix.unc["SB0571",1:2]<-c(0.5,0.5)
    oc.m.matrix.unc["GB103",1:2]<-c(0.5,0.5)
    oc.m.matrix.unc["GB153",1:2]<-c(0.5,0.5)

# FEMALE OCELLI
    oc.f.matrix.unc <- to.matrix(oc.f[tree$tip.label],levels(oc.f))
    oc.f.matrix.unc
    
    oc.f.matrix.unc["SB0280",1:2]<-c(0.5,0.5)
    oc.f.matrix.unc["SB0470",1:2]<-c(0.5,0.5)
    oc.f.matrix.unc["SB0476",1:2]<-c(0.5,0.5)
    oc.f.matrix.unc["SB0480",1:2]<-c(0.5,0.5)
    oc.f.matrix.unc["SB0529",1:2]<-c(0.5,0.5)
    oc.f.matrix.unc["SB0530",1:2]<-c(0.5,0.5)
    oc.f.matrix.unc["SB0559",1:2]<-c(0.5,0.5)
    oc.f.matrix.unc["GB111",1:2]<-c(0.5,0.5)

    
########################################################
    
### Modeltest 
    

# fitDiscrete
test.ER <- fitDiscrete(tree, wing.m, model ="ER")
plot(test.ER)
test.ARD <- fitDiscrete(tree, wing.m, model ="ARD")
plot(test.ARD)
test.SYM <- fitDiscrete(tree, wing.m, model ="SYM")
plot(test.SYM)

# irreversible model transition matrix (IRR model)
myQ<-matrix(c(0,0,1,1,0,1,1,0,0), nrow = 3,ncol=3)
rownames(myQ)<-colnames(myQ)<-c("fully","no","partially")
test.IRR <- fitDiscrete(tree, wing.m, model=myQ)

# Compare fitDiscrete models
aicc<-setNames(c(test.ER$opt$aicc,test.ARD$opt$aicc,test.SYM$opt$aicc,test.IRR$opt$aicc),c("ER","ARD","SYM","IRR"))
aic.w(aicc)



# fitMk

fitMK.ER <-fitMk(tree, wing.m, model="ER", type="discrete")
fitMK.SYM <-fitMk(tree, wing.m, model="SYM", type="discrete")
fitMK.ARD <-fitMk(tree, wing.m, model="ARD", type="discrete")

# irreversible (IRR) model transition matrix
fitMK.IRR <- fitMk(tree, wing.m, model =myQ)

# IRR model with root forced to be winged
fitMK.IRR.rootwinged <- fitMk(tree, wing.m, model =myQ,pi = c(0.5,0,0.5))

# ARD model with root forced to be wingless
fitMK.ARD.rootwingless <- fitMk(tree, wing.m, model ="ARD",pi = c(0,1,0))

# compare models (although not recommended because based on different priors!!!)
fitMK.ER.aic <- AIC(fitMK.ER)
fitMK.ARD.aic <- AIC(fitMK.ARD)
fitMK.SYM.aic <- AIC(fitMK.SYM)
fitMK.IRR.aic <- AIC(fitMK.IRR)
fitMK.IRR.rootwinged.aic <- AIC(fitMK.IRR.rootwinged)
fitMK.ARD.rootwingless.aic <- AIC(fitMK.ARD.rootwingless)

aicc2<-setNames(
  c(fitMK.ER.aic,fitMK.ARD.aic,fitMK.SYM.aic,fitMK.IRR.aic,fitMK.IRR.rootwinged.aic,fitMK.ARD.rootwingless.aic),
  c("ER","ARD","SYM","IRR","IRR.w","ARD.wl"))
aic.w(aicc2)




###################################################################################

# stochastic character mapping

###################################################################################


# run the mcmc chains
trees.wings.m <- make.simmap(tree,wing.m.matrix.unc,model="ARD",Q="mcmc",nsim=300)
trees.oc.m <- make.simmap(tree,oc.m.matrix.unc,model="ARD",Q="mcmc",nsim=300)

# binary dataset (datawings)
binary.wing.m<-setNames(wing.table.m[[3]],rownames(wing.table.m))
binary.wing.m <- as.factor(binary.wing.m) # on linux system
trees.wings.m.binary <- make.simmap(tree,binary.wing.m,model="ARD",Q="mcmc",nsim=300)

# binary IRR
mybinQ <- matrix(c(0,1,0,0), nrow=2,ncol=2)  # create IRR matrix
rownames(mybinQ)<-colnames(mybinQ)<-c("wingless","winged")
trees.wings.m.binary.IRR <- make.simmap(tree,binary.wing.m,model=mybinQ,Q="mcmc",nsim=300)



# get the states at ancestral nodes
ancst.wings.m<-describe.simmap(trees.wings.m,plot=FALSE)
ancst.oc.m<-describe.simmap(trees.oc.m,plot=FALSE)
# binary
ancst.wings.m.binary<-describe.simmap(trees.wings.m.binary,plot=FALSE)
ancst.wings.m.binary.IR<-describe.simmap(trees.wings.m.binary.IRR,plot=FALSE)



###################################################################################################

# Plot trees

cols.wings<-setNames(c("darkolivegreen1","gray8","olivedrab"),levels(wing.m))
cols.oc<-setNames(c("gray","yellow"),levels(oc.f))


# males - wings
pdffn = "ARD.males.wings.pdf"
pdf(pdffn, width=40, height=120)
plot(tree)
axisPhylo(side=1)
nodelabels(pie = ancst.wings.m$ace, piecol = cols.wings, cex = 0.2)
tiplabels(pie=wing.m.matrix.unc,piecol=cols.wings,cex=0.15)
dev.off() 

# males - ocelli
pdffn = "ARD.males.ocelli.pdf"
pdf(pdffn, width=40, height=120)
plot(tree)
axisPhylo(side=1)
nodelabels(pie = ancst.oc.m$ace, piecol = cols.oc, cex = 0.2)
tiplabels(pie=oc.m.matrix.unc,piecol=cols.oc,cex=0.15)
dev.off() 




##############################################################
##############################################################


# Transitions - Boxplots

# For binary datasets ARD and IRR (datawings)

# Libraries
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)


ard <-summary(trees.wings.m.binary)
ard
irr <-summary(trees.wings.m.binary.IRR)
irr

# all transitions for ARD
ard.all<-ard$count[,1]
mean(ard.all)
max(ard.all)
min(ard.all)
quantile(ard.all)
# transitions wingless -> winged
ard.ab<- ard$count[,2]
# transitions winged -> wingless
ard.ba<- ard$count[,3]

# same for IRR
irr.all<-irr$count[,1]
irr.ab<-irr$count[,2]   # this is zero
irr.ba<-irr$count[,3]   # this is identical to irr.all

dataset <- data.frame(
  name=c(rep("ARD a->b",300), rep("ARD b->a",300),rep("ARD combined",300),rep("IRR b=>a",300)),
  value=c(ard.ab,ard.ba,ard.all,irr.ba))

ggplot(dataset, aes(x=name, y=value, fill=name)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")

t.test(ard.ges,irr.ges,paired = TRUE)
