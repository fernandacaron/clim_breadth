rm(list=ls())

library(phytools)
library(geiger)
source("https://raw.githubusercontent.com/mgharvey/ES-sim/master/R/essim.R")

#Code for test trait-dependent diversification using tip rate correlation (TRC)

#modifying the function essim to try to speed up
essim2<-function(phy, trait, nsim = 1000, es, return.es=FALSE) {
	
    if (missing(es)) { # If inverse equal splits statistics not provided, calculate it
        #rootnode <- length(phy$tip.label) + 1
        #es <- numeric(length(phy$tip.label))
        #for (i in 1:length(es)){
        #	node <- i
        #	index <- 1
        #	qx <- 0
        #	while (node != rootnode){
        #		el <- phy$edge.length[phy$edge[,2] == node]
        #		node <- phy$edge[,1][phy$edge[,2] == node]			
        #		qx <- qx + el* (1 / 2^(index-1))			
        #		index <- index + 1
        #	}
        #	es[i] <- 1/qx
        #}		
        #names(es) <- phy$tip.label
        es <- compute_es(phy);
    }
    
    # check if all tips have traits. will work whether newly computed or passed in
    if (length(trait) < length(phy$tip.label)) {
        idx <- which(!phy$tip.label %in% names(trait));
        phy <- drop.tip(phy, idx);
    }
	
	es <- log(es[phy$tip.label]) # log transform
	trait <- trait[phy$tip.label]
	
	# Pearson's correlation between log inverse equal splits statistic and trait
	res <- cor.test(es, trait, method="pearson")

	# Fit Brownian motion model to get diffusion rate and root state estimates
	# vv <- vcv.phylo(as.phylo(phy))
	# onev <- matrix(rep(1, length(trait)), nrow=length(trait), ncol=1)
	# root <- as.vector(solve(t(onev)%*% solve(vv) %*% onev) %*% (t(onev)%*% solve(vv) %*% trait))
	# rate <- as.vector((t(trait-root) %*% solve(vv) %*% (trait-root))/length(trait))
	
	# # Brownian simulations 
	# sims <- t(rmvnorm(nsim, sigma=rate*vv))
	# rownames(sims) <- rownames(vv)
	sims<-fastBM(phy, nsim=nsim)
		
	# Pearson's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(es[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
	
	# Calculate the two-tailed p value
	corr <- res$estimate
	upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
	lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	if(missing(return.es)) { # output just rho and p value
		result <- as.vector(c(corr, pval))
		names(result) <- c("rho", "P Value")
		return(result)
	} else { # output rho, p value, and list of es values
		result <- as.vector(c(corr, pval, list(es)))
		names(result) <- c("rho", "P Value", "es")
		return(result)		
	}
}

dat<-read.csv("bioclim_data12Aug20prevalences.csv")

trAm<-read.nexus("amphibia27JUL2020.nex")
trSq<-read.nexus("squamata27JUL20.nex")
trMa<-read.nexus("mammalia_node_dated_27JUL2020.nex")
trAv<-read.nexus("aves_Ericson_27_JUL_20.nex")

#calculting es statistic with the complete trees
esAm<-esSq<-esMa<-esAv<-list()
for(i in 1:length(trAm)){
	esAm[[i]]<-compute_es(trAm[[i]])
	esSq[[i]]<-compute_es(trSq[[i]])
	esMa[[i]]<-compute_es(trMa[[i]])
	esAv[[i]]<-compute_es(trAv[[i]])
}

#separating datasets
datAm<-dat[dat$spp %in% names(esAm[[1]]),]
datSq<-dat[dat$spp %in% names(esSq[[1]]),]
datMa<-dat[dat$spp %in% names(esMa[[1]]),]
datAv<-dat[dat$spp %in% names(esAv[[1]]),]

#bio1
bio01Am<-datAm$bio1
names(bio01Am)<-datAm$spp

bio01Sq<-datSq$bio1
names(bio01Sq)<-datSq$spp

bio01Ma<-datMa$bio1
names(bio01Ma)<-datMa$spp

bio01Av<-datAv$bio1
names(bio01Av)<-datAv$spp

#bio12
bio12Am<-datAm$bio12
names(bio12Am)<-datAm$spp

bio12Sq<-datSq$bio12
names(bio12Sq)<-datSq$spp

bio12Ma<-datMa$bio12
names(bio12Ma)<-datMa$spp

bio12Av<-datAv$bio12
names(bio12Av)<-datAv$spp

#niche breadths temp
nb_tempAm<-datAm$bio5-datAm$bio6
names(nb_tempAm)<-datAm$spp

nb_tempSq<-datSq$bio5-datSq$bio6
names(nb_tempSq)<-datSq$spp

nb_tempMa<-datMa$bio5-datMa$bio6
names(nb_tempMa)<-datMa$spp

nb_tempAv<-datAv$bio5-datAv$bio6
names(nb_tempAv)<-datAv$spp

#niche breadths prec
nb_precAm<-datAm$bio16-datAm$bio17
names(nb_precAm)<-datAm$spp

nb_precSq<-datSq$bio16-datSq$bio17
names(nb_precSq)<-datSq$spp

nb_precMa<-datMa$bio16-datMa$bio17
names(nb_precMa)<-datMa$spp

nb_precAv<-datAv$bio16-datAv$bio17
names(nb_precAv)<-datAv$spp

#testing the relationship between es and temperature (bio1 - annual mean temperature) and 
#precipitation (bio12 - annual precipitation) for 1000 trees with 1000 simulations for each tree
essim01_Am<-essim12_Am<-list()
for(i in 1:1000){
	essim01_Am[[i]]<-essim2(phy=trAm[[i]], trait=bio01Am, es=esAm[[i]], nsim=1000)
	essim12_Am[[i]]<-essim2(phy=trAm[[i]], trait=bio12Am, es=esAm[[i]], nsim=1000)
}

essim01_Sq<-essim12_Sq<-list()
for(i in 1:1000{
	essim01_Sq[[i]]<-essim2(phy=trSq[[i]], trait=bio01Sq, es=esSq[[i]], nsim=1000)
	essim12_Sq[[i]]<-essim2(phy=trSq[[i]], trait=bio12Sq, es=esSq[[i]], nsim=1000)
}

essim01_Ma<-essim12_Ma<-list()
for(i in 1:1000){
	essim01_Ma[[i]]<-essim2(phy=trMa[[i]], trait=bio01Ma, es=esMa[[i]], nsim=1000)
	essim12_Ma[[i]]<-essim2(phy=trMa[[i]], trait=bio12Ma, es=esMa[[i]], nsim=1000)
}

essim01_Av<-essim12_Av<-list()
for(i in 1:1000){
	essim01_Av[[i]]<-essim2(phy=trAv[[i]], trait=bio01Av, es=esAv[[i]], nsim=1000)
	essim12_Av[[i]]<-essim2(phy=trAv[[i]], trait=bio12Av, es=esAv[[i]], nsim=1000)
}

#testing the relationship between es and niche breadths temperature (bio5 - bio6) and 
#precipitation (bio16 - bio17) for 1000 trees with 1000 simulations for each tree
essim_nb_tempAm<-essim_nb_precAm<-list()
for(i in 1:1000){
	essim_nb_tempAm[[i]]<-essim2(phy=trAm[[i]], trait=nb_tempAm, es=esAm[[i]], nsim=1000)
	essim_nb_precAm[[i]]<-essim2(phy=trAm[[i]], trait=nb_precAm, es=esAm[[i]], nsim=1000)
}

essim_nb_tempSq<-essim_nb_precSq<-list()
for(i in 1:1000){
	essim_nb_tempSq[[i]]<-essim2(phy=trSq[[i]], trait=nb_tempSq, es=esSq[[i]], nsim=1000)
	essim_nb_precSq[[i]]<-essim2(phy=trSq[[i]], trait=nb_precSq, es=esSq[[i]], nsim=1000)
}

essim_nb_tempMa<-essim_nb_precMa<-list()
for(i in 1:1000){
	essim_nb_tempMa[[i]]<-essim2(phy=trMa[[i]], trait=nb_tempMa, es=esMa[[i]], nsim=1000)
	essim_nb_precMa[[i]]<-essim2(phy=trMa[[i]], trait=nb_precMa, es=esMa[[i]], nsim=1000)
}

essim_nb_tempAv<-essim_nb_precAv<-list()
for(i in 1:1000){
	essim_nb_tempAv[[i]]<-essim2(phy=trAv[[i]], trait=nb_tempAv, es=esAv[[i]], nsim=1000)
	essim_nb_precAv[[i]]<-essim2(phy=trAv[[i]], trait=nb_precAv, es=esAv[[i]], nsim=1000)
}

essim01_Am<-do.call(rbind, essim01_Am)
essim01_Sq<-do.call(rbind, essim01_Sq)
essim01_Ma<-do.call(rbind, essim01_Ma)
essim01_Av<-do.call(rbind, essim01_Av)

essim12_Am<-do.call(rbind, essim12_Am)
essim12_Sq<-do.call(rbind, essim12_Sq)
essim12_Ma<-do.call(rbind, essim12_Ma)
essim12_Av<-do.call(rbind, essim12_Av)

essim_nb_tempAm<-do.call(rbind, essim_nb_tempAm)
essim_nb_tempSq<-do.call(rbind, essim_nb_tempSq)
essim_nb_tempMa<-do.call(rbind, essim_nb_tempMa)
essim_nb_tempAv<-do.call(rbind, essim_nb_tempAv)

essim_nb_precAm<-do.call(rbind, essim_nb_precAm)
essim_nb_precSq<-do.call(rbind, essim_nb_precSq)
essim_nb_precMa<-do.call(rbind, essim_nb_precMa)
essim_nb_precAv<-do.call(rbind, essim_nb_precAv)

pdf("FigureS2.pdf")
m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = TRUE)

layout(mat = m,heights = c(0.3,0.3,0.1))

hist(essim01_Am[,1], xlab="rho", main="BIO1", xlim=c(-0.22, 0.13), ylim=c(0,260), breaks=seq(-0.22, 0.13, 0.01), col=rgb(92/255,18/255,110/255, 0.5), border=inferno(20)[6])
hist(essim01_Sq[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(177/255, 50/255, 90/255, 0.5), border=inferno(20)[10])
hist(essim01_Ma[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(241/255, 112/255, 32/255, 0.5), border=inferno(20)[14])
hist(essim01_Av[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(254/255, 183/255, 45/255, 0.4), border=plasma(20)[17])

hist(essim12_Am[,1], xlab="rho", main="BIO12", xlim=c(-0.22, 0.13), ylim=c(0,270), breaks=seq(-0.22, 0.13, 0.01), col=rgb(92/255,18/255,110/255, 0.5), border=inferno(20)[6])
hist(essim12_Sq[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(177/255, 50/255, 90/255, 0.5), border=inferno(20)[10])
hist(essim12_Ma[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(241/255, 112/255, 32/255, 0.5), border=inferno(20)[14])
hist(essim12_Av[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(254/255, 183/255, 45/255, 0.4), border=plasma(20)[17])

hist(essim_nb_tempAm[,1], xlab="rho", main="Temperature niche breadth", ylim=c(0,250), xlim=c(-0.22, 0.13), breaks=seq(-0.22, 0.13, 0.01), col=rgb(92/255,18/255,110/255, 0.5), border=inferno(20)[6])
hist(essim_nb_tempSq[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(177/255, 50/255, 90/255, 0.5), border=inferno(20)[10])
hist(essim_nb_tempMa[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(241/255, 112/255, 32/255, 0.5), border=inferno(20)[14])
hist(essim_nb_tempAv[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(254/255, 183/255, 45/255, 0.4), border=plasma(20)[17])

hist(essim_nb_precAm[,1], xlab="rho", main="Precipitation niche breadth", xlim=c(-0.22, 0.13), ylim=c(0,300), breaks=seq(-0.22, 0.13, 0.01), col=rgb(92/255,18/255,110/255, 0.5), border=inferno(20)[6])
hist(essim_nb_precSq[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(177/255, 50/255, 90/255, 0.5), border=inferno(20)[10])
hist(essim_nb_precMa[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(241/255, 112/255, 32/255, 0.5), border=inferno(20)[14])
hist(essim_nb_precAv[,1], breaks=seq(-0.22, 0.13, 0.01), add=T, col=rgb(254/255, 183/255, 45/255, 0.4), border=plasma(20)[17])

par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=6,legend = c("Amphibia", "Squamata", "Mammals", "Birds"), col=c(rgb(92/255,18/255,110/255, 0.7), rgb(177/255, 50/255, 90/255, 0.7),
		rgb(241/255, 112/255, 32/255, 0.7), rgb(254/255, 183/255, 45/255, 0.6)), border=c(inferno(20)[6],inferno(20)[10],inferno(20)[14], plasma(20)[17]), 
		bty = 'n', pch=15, cex=1.2, text.font=3)
		
dev.off()

#organizing the data to put in tables
essim01_Am<-cbind(essim01_Am, rep("Amphibia",1000))
essim01_Sq<-cbind(essim01_Sq, rep("Squamata",1000))
essim01_Ma<-cbind(essim01_Ma, rep("Mammals",1000))
essim01_Av<-cbind(essim01_Av, rep("Birds",1000))

essim12_Am<-cbind(essim12_Am, rep("Amphibia",1000))
essim12_Sq<-cbind(essim12_Sq, rep("Squamata",1000))
essim12_Ma<-cbind(essim12_Ma, rep("Mammals",1000))
essim12_Av<-cbind(essim12_Av, rep("Birds",1000))

essim_nb_tempAm<-cbind(essim_nb_tempAm, rep("Amphibia",1000))
essim_nb_tempSq<-cbind(essim_nb_tempSq, rep("Squamata",1000))
essim_nb_tempMa<-cbind(essim_nb_tempMa, rep("Mammals",1000))
essim_nb_tempAv<-cbind(essim_nb_tempAv, rep("Birds",1000))

essim_nb_precAm<-cbind(essim_nb_precAm, rep("Amphibia",1000))
essim_nb_precSq<-cbind(essim_nb_precSq, rep("Squamata",1000))
essim_nb_precMa<-cbind(essim_nb_precMa, rep("Mammals",1000))
essim_nb_precAv<-cbind(essim_nb_precAv, rep("Birds",1000))

essim_01<-as.data.frame(rbind(essim01_Am, essim01_Sq, essim01_Ma, essim01_Av), stringsAsFactors=FALSE)
essim_12<-as.data.frame(rbind(essim12_Am, essim12_Sq, essim12_Ma, essim12_Av), stringsAsFactors=FALSE)
essim_nb_temp<-as.data.frame(rbind(essim_nb_tempAm, essim_nb_tempSq, essim_nb_tempMa, essim_nb_tempAv), stringsAsFactors=FALSE)
essim_nb_prec<-as.data.frame(rbind(essim_nb_precAm, essim_nb_precSq, essim_nb_precMa, essim_nb_precAv), stringsAsFactors=FALSE)

essim<-cbind(essim_01, essim_12, essim_nb_temp, essim_nb_prec)
essim<-as.data.frame(as.matrix(essim), stringsAsFactors=FALSE)
colnames(essim)<-c("rho_01", "p_01", "taxon_01","rho_12", "p_12", "taxon_12",
					"rho_nb_temp", "p_nb_temp", "taxon_nb_temp", "rho_nb_prec", 
					"p_nb_prec", "taxon_nb_prec")


#getting the 95% confidence intervals from the estimates
l95<-function(x) quantile(x, probs=0.025)
h95<-function(x) quantile(x, probs=0.975)
​
final<-data.frame(
median_rho_BIO1=aggregate(as.numeric(essim$rho_01), list(essim$taxon_01), median),
low95_rho_BIO1=aggregate(as.numeric(essim$rho_01), list(essim$taxon_01), l95)[,2],
high95_rho_BIO1=aggregate(as.numeric(essim$rho_01), list(essim$taxon_01), h95)[,2],

median_p_BIO1=aggregate(as.numeric(essim$p_01), list(essim$taxon_01), median)[,2],
low95_p_BIO1=aggregate(as.numeric(essim$p_01), list(essim$taxon_01), l95)[,2],
high95_p_BIO1=aggregate(as.numeric(essim$p_01), list(essim$taxon_01), h95)[,2],

median_rho_BIO12=aggregate(as.numeric(essim$rho_12), list(essim$taxon_12), median)[,2],
low95_rho_BIO12=aggregate(as.numeric(essim$rho_12), list(essim$taxon_12), l95)[,2],
high95_rho_BIO12=aggregate(as.numeric(essim$rho_12), list(essim$taxon_12), h95)[,2],

median_p_BIO12=aggregate(as.numeric(essim$p_12), list(essim$taxon_12), median)[,2],
low95_p_BIO12=aggregate(as.numeric(essim$p_12), list(essim$taxon_12), l95)[,2],
high95_p_BIO12=aggregate(as.numeric(essim$p_12), list(essim$taxon_12), h95)[,2],

median_rho_NBTEMP=aggregate(as.numeric(essim$rho_nb_temp), list(essim$taxon_nb_temp), median)[,2],
low95_rho_NBTEMP=aggregate(as.numeric(essim$rho_nb_temp), list(essim$taxon_nb_temp), l95)[,2],
high95_rho_NBTEMP=aggregate(as.numeric(essim$rho_nb_temp), list(essim$taxon_nb_temp), h95)[,2],

median_p_NBTEMP=aggregate(as.numeric(essim$p_nb_temp), list(essim$taxon_nb_temp), median)[,2],
low95_p_NBTEMP=aggregate(as.numeric(essim$p_nb_temp), list(essim$taxon_nb_temp), l95)[,2],
high95_p_NBTEMP=aggregate(as.numeric(essim$p_nb_temp), list(essim$taxon_nb_temp), h95)[,2],

median_rho_NBPREC=aggregate(as.numeric(essim$rho_nb_prec), list(essim$taxon_nb_prec), median)[,2],
low95_rho_NBPREC=aggregate(as.numeric(essim$rho_nb_prec), list(essim$taxon_nb_prec), l95)[,2],
high95_rho_NCPREB=aggregate(as.numeric(essim$rho_nb_prec), list(essim$taxon_nb_prec), h95)[,2],

median_p_NBPREC=aggregate(as.numeric(essim$p_nb_prec), list(essim$taxon_nb_prec), median)[,2],
low95_p_NBPREC=aggregate(as.numeric(essim$p_nb_prec), list(essim$taxon_nb_prec), l95)[,2],
high95_p_NCPREB=aggregate(as.numeric(essim$p_nb_prec), list(essim$taxon_nb_prec), h95)[,2]
)

rownames(final)<-final[,1]
final<-final[c("Amphibia", "Squamata", "Mammalia", "Aves"),-1]
​
write.csv(final, "Table_s5_unformatted.csv")
