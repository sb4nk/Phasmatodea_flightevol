########################################################################################

# R code used in Bank & Bradler (2022)
# https://github.com/sb4nk/Phasmatodea_flightevol

# Estimation of the phylogenetic signal using D statistics and Pagel's lambda, and 
# the comparison of the observed number of evolutionary transitions for each trait 
# against a randomised character matrix.

########################################################################################


library(ape)
library(phytools)
library(geiger)
library(picante)
library(caper)
library("phangorn")


tree <- read.nexus("treefile.BI_constraintsB2.nex")


### CAPER #### calculate D statistic

### WINGS
wing.table.caper <- read.csv("morphdata_males.txt")[ ,c(1,4)]
wing.caper <- comparative.data(phy=tree, data=wing.table.caper, names.col = ID)
res.wing <- phylo.d(data=wing.caper, binvar = datawings, permut= 10000)
res.wing

#### FLIGHT
flight.table.caper <- read.csv("morphdata_males.txt")[ ,c(1,5)]
flight.caper <- comparative.data(phy=tree, data=flight.table.caper, names.col = ID)
res.flight <- phylo.d(data=flight.caper, binvar = dataflight, permut= 10000)
res.flight

### OCELLI
ocelli.table.caper <- read.csv("morphdata_males.txt")[ ,c(1,3)]
ocelli.caper <- comparative.data(phy=tree, data=ocelli.table.caper, names.col = ID)
res.ocelli <- phylo.d(data=ocelli.caper, binvar = ocelli, permut= 10000)
res.ocelli



#############################################################################################################
#############################################################################################################



# Pagel's lambda


table <- read.csv("morphdata_males.txt", row.names = 1)


# Is there strong phylo singal in presence of wings? ############################################

wings <- table$datawings
names(wings) <- row.names(table)
wings <- as.factor(wings)

# exchange 0 for A (wings absent)and 1 for B (wings present)
table$datawings[table$datawings == "0"] <- "A"   
table$datawings[table$datawings == "1"] <- "B"

wings.renamed <- table$datawings
names(wings.renamed) <- row.names(table)
wings.renamed <- as.factor(wings.renamed)

# Is there strong phylo singal in presence of flight (fully-winged species)? ####################

flight <- table$dataflight
names(flight) <- row.names(table)
flight <- as.factor(flight)

# exchange 0 for A (flightless)and 1 for B (flighted)
table$dataflight[table$dataflight == "0"] <- "A"   
table$dataflight[table$dataflight == "1"] <- "B"

flight.renamed <- table$dataflight
names(flight.renamed) <- row.names(table)
flight.renamed <- as.factor(flight.renamed)


# Is there strong phylo singal in presence ocelli ?####################
ocelli <- table$ocelli
names(ocelli) <- row.names(table)
ocelli <- as.factor(ocelli)



# Using fitDiscrete. Pagel's lambda is a multiplier of the off-diagonal elements 
# of a variance-covariance matrix, which best fits the distribution of data
# at the tips of the phylogeny. Values vary between 0 (phylogenetic independence)
# and 1 (traits evolve according to a Brownian process)
wings_lambda <- fitDiscrete(tree, wings.renamed,transform="lambda", model="ARD")
wings_lambda
wings_lambda$opt$lambda

# flight
flight_lambda <- fitDiscrete(tree, flight.renamed, transform="lambda", model="ARD")
flight_lambda
flight_lambda$opt$lambda

# ocelli
ocelli_lambda <- fitDiscrete(tree, ocelli, transform="lambda", model="ARD")
ocelli_lambda
ocelli_lambda$opt$lambda




########################################################################################

# Randomizations of a character on a tree

# from https://github.com/juliema/publications/blob/master/BrueeliaMS/Maddison.Slatkin.R
# Will tell you whether there are LESS or MORE evolutionary transitions than expected by 
# randomization process. 

phylo.signal.disc <- function(trait,phy, rep = 999,cost=NULL){
  lev <- attributes(factor(trait))$levels
  if (length(lev) == length(trait))
    stop("Are you sure this variable is categorical?")
  if(is.null(cost)){
    cost1 <- 1-diag(length(lev))
  }
  else {
    if (length(lev) != dim(cost)[1])
      stop("Dimensions of the character state transition matrix do not agree with the number of levels")
    cost1<- t(cost)
  }
  dimnames(cost1) <- list(lev,lev)
  trait <- as.numeric(trait)
  attributes(trait)$names <- phy$tip
  NULL.MODEL <- matrix(NA,rep,1)
  obs <- t(data.frame(trait))
  obs <- phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
  OBS <- parsimony(phy,obs,method="sankoff",cost=cost1)
  for (i in 1:rep){
    null <- sample(as.numeric(trait))
    attributes(null)$names <- attributes(trait)$names
    null <- t(data.frame(null))
    null <- phyDat(t(null),type="USER",levels=attributes(factor(null))$levels)
    NULL.MODEL[i,]<-parsimony(phy,null,method="sankoff",cost=cost1)
    P.value <- sum(OBS >= NULL.MODEL)/(rep + 1)
    #print(OBS)
    #print(NULL.MODEL)
  }
  par(mfrow=c(1,2))
  hist(NULL.MODEL,xlab="Transitions.in.Randomizations",xlim=c(min(c(min(NULL.MODEL,OBS-1))),max(NULL.MODEL)+1))
  arrows(OBS,rep/10,OBS,0,angle=20,col="red",lwd=4)
  phy$tip.label <- rep(".",length(trait))
  plot(phy,tip.col=trait+10,cex=250/length(trait),font=1)
  title("Character states")
  par(mfrow=c(1,1))
  
  OUTPUT1 <- t(data.frame(Number.of.Levels = length(attributes(factor(trait))$levels), Evolutionary.Transitions.Observed=OBS,Evolutionary.Transitions.Randomization.Median=median(NULL.MODEL),Evolutionary.Transitions.Randomization.Min=min(NULL.MODEL),Evolutionary.Transitions.Randomization.Max=max(NULL.MODEL),P.value))
  
  if(is.null(cost)){
    list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.UNORDERED.PARSIMONY = t(cost1), .DATA=NULL.MODEL, .OBSERVED=OBS)
  }
  else {
    list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.FROM.ROW.TO.COL = t(cost1))        }
}


# wings
wing <- phylo.signal.disc(wings.renamed, tree)
wing$.Randomization.Results 
#wing$.DATA[,1]
#wing$.OBSERVED
mean(wing$.DATA[,1])


# flight
flight <- phylo.signal.disc(flight.renamed, tree)
flight$.Randomization.Results
mean(flight$.DATA[,1])


# ocelli
ocelli <- phylo.signal.disc(ocelli, tree)
ocelli$.Randomization.Results 
mean(ocelli$.DATA[,1])
