#setup packages

install.packages("ape")
install.packages("phytools")

library(ape)
library(phytools)


#load in the phylogeny
tree <- read.tree("RAxML_bipartitions.roary_for_europe_final_2019.tre")
plot(tree)

#reroot the tree by the outgroup
roottree <- root(tree, outgroup = "IVIA5235_SP")

#lets check out what our tree file contains
head(tree$tip.label)
head(tree$edge.length) #tab complete with an object and $ to see possibilities

#making a new variable based on the taxa included that we will use
Species <- tree$tip.label
write.csv(Species,file = "xf_Species.csv")
xf_data <- read.csv("xf_data.csv", row.names = 1)
char_data <- as.matrix(xf_data)[,1]
char_data <-setNames(xf_data[,1],rownames(xf_data))

#lets see what our rerooted tree looks like with a different plot function
plotTree(roottree, fsize=0.3, lwd=.8, offset = .5) 



#ancestral state reconstructions

#run the equal rates model
fitER <- rerootingMethod(roottree, char_data, model="ER")

#Check out our equal rates model
head(fitER$marginal.anc)

#run the symmetrical rates model
fitSYM <- rerootingMethod(roottree, char_data, model="SYM")

#Check out our sym rates model
head(fitSYM$marginal.anc)

#lets look under the hood
print("Equal rates matrix")
fitER$Q

print("Symmetrical rates matrix")
fitSYM$Q

#plotting ER

#assign a color set
cols<-rainbow(n = 6)

#plot tree
plotTree(roottree, fsize=0.3, lwd=.8, offset = .5) 

#add labels to the nodes that come from our reconstruction results file
nodelabels(pie = fitER$marginal.anc, piecol=cols, cex=0.3) 
legend(x=0, y=10, legend = c("BR", "FR", "IT","SP", "US_CA", "US_SE"), lty=1, col=cols, cex=.8)

#log likelihood
fitER$loglik

#AIC score
AICer <- -2*fitER$loglik + 2*1
AICer




#plotting sym

#assign a color set
cols<-rainbow(n = 6)

#plot tree
plotTree(roottree, fsize=0.3, lwd=.8, offset = .5) 

#add labels to the nodes that come from our reconstruction results file
nodelabels(pie = fitSYM$marginal.anc, piecol=cols, cex=0.3) 
legend(x=0, y=10, legend = c("BR", "FR", "IT","SP", "US_CA", "US_SE"), lty=1, col=cols, cex=.8)

#log likelihood
fitSYM$loglik

#AIC score
AICsym <- -2*fitSYM$loglik + 2*6
AICsym

