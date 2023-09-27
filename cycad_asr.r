## load package and time calibrated phylogeny
library(phytools)
phy <- read.nexus("C:\\Users\\whbri\\Desktop\\Research\\Cycad\\Cycadales_TED_divFBD_IGR_fullconstraints.con.tre")

## compile n-fixing dataset
root <- max(nodeHeights(phy),root.edge=TRUE)
Nfix <- data.frame("tip" = phy$tip.label,
					"age.bp" = round(sapply(1:Ntip(phy),function(x) nodeheight(phy,x,root.edge=FALSE)-root),1)
					)
Nfix$Y <- 0.5
Nfix[Nfix$age.bp == 0,]$Y <- 1	## designate all extant taxa as N fixing

## compile fossil n-fixing data
Nfix.fossil <- read.csv("C:\\Users\\whbri\\Desktop\\Research\\Cycad\\fossil-N.csv")
nitroF <- Nfix.fossil[Nfix.fossil$Y == 1,]$tip
normlF <- Nfix.fossil[Nfix.fossil$Y == 0,]$tip

Nfix[Nfix$tip %in% nitroF,]$Y <- 1
Nfix[Nfix$tip %in% normlF,]$Y <- 0

Nfix$N <- 1-Nfix$Y

P <- as.matrix(Nfix[,c("Y","N")])
rownames(P) <- Nfix$tip

## run ancestral state estimates
er <- make.simmap(phy,P,model = "ER", nsim = 1)
ard <- make.simmap(phy,P,model = "ARD", nsim = 1)
AIC(er,ard)

model <- make.simmap(phy,P,Q = ard$Q, nsim = 1000)
smry <- summary(model)


	## plot reconstruction without tip pies for uninformative taxa
		#pdf(height = 18, width = 9,file.choose())
	plot(phy,cex = 0.45,label.offset = 3)
	inf <- which(!rownames(smry$tips) %in% Nfix[Nfix$Y == 0.5,]$tip)
	nodelabels(pie=smry$ace,piecol=c("black","white"),cex=0.20)
	tiplabels(tip = inf,pie=smry$tips[inf,],piecol=c("black","white"),cex=0.15)
	abline(v = mx-c(0,23,66,145,201,252,299,359))


## trim dataset/phylogeny to only have a single extant representative for clades
Nfix$Genus <- gsub("_.*","",Nfix$tip)
ext <- unique(Nfix[Nfix$age.bp == 0,]$Genus)

set.seed(12345)
ext_incl <-
lapply(ext, function(x) {
	which(Nfix$Genus == x & Nfix$age.bp == 0)[[1]]
	})


plotN <- Nfix[c(unlist(ext_incl),which(Nfix$age.bp < 0 )),]
phy$node.label <- 1:Nnode(phy)
plotPhy <-  drop.tip(phy,setdiff(phy$tip.label,plotN$tip))

	## plot reconstruction without tip pies for uninformative taxa
		#pdf(height = 13, width = 9,file.choose())
	plot(plotPhy,cex = 0.65,label.offset = 3)
	inf <- which(rownames(smry$tips) %in% plotN[!plotN$Y == 0.5,]$tip)
	nodelabels(pie=smry$ace[plotPhy$node.label,],piecol=c("black","white"),cex=0.35)
	tiplabels(tip = which(plotPhy$tip.label %in% rownames(smry$tips[inf,])),pie=smry$tips[inf,],piecol=c("black","white"),cex=0.3)
	#abline(v = mx-c(0,23,66,145,201,252,299,359))
