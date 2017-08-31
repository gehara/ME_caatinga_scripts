library(gbm)
library(caret)
library(pROC)
library(abc)
setwd("~/Github/Networksimulator/R")
source("network.sim.R")

######### observed distance
setwd("/home/marcelo/Dropbox/AMNH/Lampropeltis/gene_fasta_11taxa/")
beast.tree<-"(((abnorma:2.413413791,micropholis:2.413413791):2.046516169,polyzona:4.459929961):2.542542814,((alterna:4.554400812,((zonata:1.534640091,multifasciata:1.534640091):1.423238548,(knoblochi:1.634034342,pyromelana:1.634034342):1.323844296):1.596522174):1.244838359,(webbi:2.892264576,(mexicana:1.539058445,ruthveni:1.539058445):1.353206132):2.906974595):1.203233603);"

observed.dist<-get.observed(tree.b=get.tree.info(beast.tree),path=getwd(),spectral=F, distance=T, moments=T)
observed.dist<-observed.dist[[1]]

### get informatin about tree tips
get.tree.info(beast.tree)
### beast plot
par(mfrow=c(1,2))
plot(read.tree(text=get.tree.info(beast.tree)[[2]][1]))
axisPhylo()
plot(read.tree(text=beast.tree),direction = "leftwards")
axisPhylo()

### beast.priors

# "9 11"  "2.892264576"
#  "2 3"   "4.45992996" 

hib.priors<-c(2.9,4.4,4.4,0,0.9) # time of split: LOWER,UPPER; minor node UPPER; minor porportion LOWER,UPPER;
hib.clade<-c(9,10,11)
major.sister<-c(4,5,6,7,8)
minor.sister<-3

# segsites: if true calculate a tree with the simualted segregating sites.
# coaltrees: use the simulated coalescent tree
# distance: get the distances as output
# spectral: get spectral eigenvalues from Panda package

setwd("~/Github/PhyNeSiCo/")

write.table(mig,"mig.txt",quote=F, row.names = F,col.names = T, sep="\t", append = F)
#write.table(bif.beast,"bif_beast.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
write.table(hyb,"hyb.txt",quote=F, row.names = F,col.names = T, sep="\t", append = F)

for(i in 1:9){
  
  #  bif<-sim.sp.tree(tree=beast.tree, Ne.prior = c(100000,2000000),time.modif=c(0.7,1.6),
  #                        nsims=1000,bp=c(1420.4,302.5),mi=c(1e-10,4e-10),gen.time=1,
  #                       biffurcating = T,ntrees=304,segsites=T,coaltrees = F,
  #                      time.scalar=1000000, distance=T,spectral=F)
  
  
  hyb<-sim.sp.tree(tree=beast.tree, Ne.prior = c(100000,2000000),time.modif=c(0.7,1.6),
                   hib.clade = hib.clade, hib.priors = hib.priors, major.sister = major.sister, 
                   minor.sister = minor.sister,
                   nsims=1000,bp=c(1420.4,302.5),mi=c(1e-10,4e-10),gen.time=1,
                   biffurcating = F, admixture=T, migration=F,
                   ntrees=304,segsites=T,coaltrees = F,
                   time.scalar=1000000, distance=T,spectral=F)
  
  mig<-sim.sp.tree(tree=beast.tree, Ne.prior = c(100000,2000000),time.modif=c(0.7,1.6),
                   hib.clade = hib.clade, hib.priors = hib.priors, major.sister = major.sister, 
                   minor.sister = minor.sister,
                   nsims=1000,bp=c(1420.4,302.5),mi=c(1e-10,4e-10),gen.time=1,mig=c(0,1),
                   biffurcating = F, admixture=F, migration=T,
                   ntrees=304,segsites=T,coaltrees = F,
                   time.scalar=1000000, distance=T,spectral=F)
  
  
  write.table(mig,"mig.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  #write.table(bif.beast,"bif_beast.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  write.table(net,"hyb.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  
}

setwd("~/Github/Networksimulator/")
mig<-read.table("mig_beast.txt",sep="\t",header=T)
hyb<-read.table("net_beast.txt",sep="\t",header=T)
bif<-read.table("bif_beast.txt",sep="\t",header=T)

colnames(bif.beast[,6:ncol(bif.beast)])->colnames(observed.dist)

goodness<-gfit(observed.dist[,1:55],bif.beast[,6:60],nb.replicate = 100,tol=1)
plot(goodness)
summary(goodness)

goodness<-gfit(observed.dist[1:55],hyb[,8:62],nb.replicate = 100,tol=1)
plot(goodness)
summary(goodness)

goodness<-gfit(observed.dist[1:55],mig[,8:62],nb.replicate = 100,tol=1)
plot(goodness)
summary(goodness)
head(mig)


index<-c(rep("bif",1000),
         rep("mig",nrow(mig)),
         rep("hyb",nrow(hyb)))

#models.astral<-rbind(bif.astral[,6:ncol(bif.astral)],net.astral[,6:ncol(net.astral)])

models<-rbind(bif[1:1000,6:60],#ncol(bif.beast)],
              mig[1:1000,8:62],#:ncol(mig.beast)],
              hyb[1:1000,8:62])#ncol(net.beast)])

models<-data.frame(models)

observed.dist<-t(data.frame(observed.dist))
colnames(models)->colnames(observed.dist)

head(models)

descrCor <- cor(models)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .90)
models.uncor <- models[,-highlyCorDescr]
head(models.uncor)


observed.dist<-observed.dist[,1:55]
head(observed.dist)

observed.uncor<-observed.dist[-highlyCorDescr]
t(data.frame(observed.uncor))->observed.uncor
colnames(models.uncor)->colnames(observed.uncor)

head(models.uncor)
head(observed.uncor)


################### all variables neuronet
##########################################
##########################################
models<-cbind(models,index)
outcomeName <- 'index'
predictorsNames <- names(models)[names(models) != outcomeName]
predictorsNames <- predictorsNames[grep("mean",predictorsNames)]
splitIndex <- createDataPartition(models[,outcomeName], p = .75, list = FALSE, times = 1)
train <- models[ splitIndex,]
test  <- models[-splitIndex,]
objControl <- trainControl(method='boot', number=10, returnResamp='final', 
                           #summaryFunction = twoClassSummary,
                           classProbs = TRUE)
nnetModel <-train(train[,predictorsNames], train[,outcomeName], 
                  method="nnet",maxit=2000, 
                  trControl=objControl,  
                  metric = "Accuracy",
                  preProc = c("center", "scale"))
predictions <- predict(object=nnetModel, test[,predictorsNames], type='raw')
accu<-postResample(pred=predictions, obs=as.factor(test[,outcomeName]))[1]

observed.dist<-t(data.frame(observed.dist[1:55]))
colnames(models[,1:55])->colnames(observed.dist)
predict(object=nnetModel, data.frame(observed.dist), type='prob')


##############################################
models.uncor<-cbind(models.uncor,index)
outcomeName <- 'index'
predictorsNames <- names(models.uncor)[names(models.uncor) != outcomeName]
splitIndex <- createDataPartition(models.uncor[,outcomeName], p = .75, list = FALSE, times = 1)
train <- models.uncor[ splitIndex,]
test  <- models.uncor[-splitIndex,]
objControl <- trainControl(method='boot', number=20, returnResamp='final', 
                           #summaryFunction = twoClassSummary, 
                           classProbs = TRUE)
nnetModel <-train(train[,predictorsNames], train[,outcomeName], 
                  method="nnet",maxit=1000, 
                  trControl=objControl,  
                  metric = "Accuracy",
                  preProc = c("center", "scale"))
predictions <- predict(object=nnetModel, test[,predictorsNames], type='raw')
accu<-postResample(pred=predictions, obs=as.factor(test[,outcomeName]))[1]

predict(object=nnetModel, data.frame(observed.uncor), type='prob')
###############################################################################
models2<-cbind(models2,index)
outcomeName <- 'index'
predictorsNames <- names(models2)[names(models2) != outcomeName]
splitIndex <- createDataPartition(models2[,outcomeName], p = .75, list = FALSE, times = 1)
train <- models2[ splitIndex,]
test  <- models2[-splitIndex,]
objControl <- trainControl(method='boot', number=20, returnResamp='final', 
                           #summaryFunction = twoClassSummary, 
                           classProbs = TRUE)
nnetModel <-train(train[,predictorsNames], train[,outcomeName], 
                  method="nnet",maxit=1000, 
                  trControl=objControl,  
                  metric = "Accuracy",
                  preProc = c("center", "scale"))
predictions <- predict(object=nnetModel, test[,predictorsNames], type='raw')
accu<-postResample(pred=predictions, obs=as.factor(test[,outcomeName]))[1]

obs2<-data.frame(obs2)
colnames(obs2)<-colnames(models2)[1:2]
predict(object=nnetModel, obs2, type='prob')

#
##############################################################################
##############################################################################
library(ggplot2)

index<-c(rep("astral",nrow(bif.astral)),rep("beast",nrow(bif.beast)))
models<-rbind(bif.astral[,2:ncol(bif.astral)],bif.beast[,2:ncol(bif.beast)])

colnames(models)->colnames(observed.spec)
xx<-rbind(models,observed.spec)

data<-c(index,rep("observed",nrow(observed.spec)))

PCA<-prcomp(xx, center = T, scale. = T, retx=T)
plot(PCA)
scores <- data.frame(PCA$x[,1:10])
ggplot(scores, aes(x=PC1, y=PC3))+
  theme(legend.text = element_text(size = 16, face = "bold"))+
  geom_point(aes(colour=data, size=data, shape=data))+
  scale_size_manual(values=c(2,2,1))+
  scale_colour_brewer(palette="Set1")
dev.off()

head(hyb)
bif2<-bif[1:1000,6:60]
hyb2<-hyb[1:1000,8:62]
mig2<-mig[1:1000,8:62]

colnames(bif2)
bif2[,c(8:10,17:19,25:27)]
bif.minor.hybrid<-rowMeans(bif2[,c(25:27)])
bif.major.hybrid<-rowMeans(bif2[,c(32:34,38:40,43:45,47:52)])
bif2<-cbind(bif.major.hybrid,bif.minor.hybrid)

hyb.minor.hybrid<-rowMeans(hyb2[,c(25:27)])
hyb.major.hybrid<-rowMeans(hyb2[,c(32:34,38:40,43:45,47:52)])
hyb2<-cbind(hyb.major.hybrid,hyb.minor.hybrid)

mig.minor.hybrid<-rowMeans(mig2[,c(25:27)])
mig.major.hybrid<-rowMeans(mig2[,c(32:34,38:40,43:45,47:52)])
mig2<-cbind(mig.major.hybrid,mig.minor.hybrid)

plot(bif.minor.hybrid~bif.major.hybrid)
points(hyb.minor.hybrid~hyb.major.hybrid,col="red")
points(mig.minor.hybrid~mig.major.hybrid,col="gray")
points(obs.minor.hybrid~obs.major.hybrid,col="yellow")

obs.minor.hybrid<-mean(observed.dist[c(25:27)])
obs.major.hybrid<-mean(observed.dist[c(32:34,38:40,43:45,47:52)])
obs2<-c(obs.major.hybrid,obs.minor.hybrid)

par(mfrow=c(3,2))
plot.sim.obs(mig2,obs2)
plot.sim.obs(hyb2,obs2)
plot.sim.obs(bif2,obs2)

hist(mig2[, 1], breaks = 20, xlab = "major-hybrid distance", main = "",xlim=c(0.002,0.01))
abline(v = obs2[1], col = 2)
hist(mig2[, 2], breaks = 10, xlab = "major-minor distance", main = "",xlim=c(0.002,0.01),ylim=c(0,5000))
abline(v = obs2[2], col = 2)

hist(hyb2[, 1], breaks = 20, xlab = "major-hybrid distance", main = "",xlim=c(0.002,0.01))
abline(v = obs2[1], col = 2)
hist(hyb2[, 2], breaks = 10, xlab = "major-minor distance", main = "",xlim=c(0.002,0.01),ylim=c(0,5000))
abline(v = obs2[2], col = 2)

hist(bif2[, 1], breaks = 20, xlab = "major-hybrid distance", main = "",xlim=c(0.002,0.01))
abline(v = obs2[1], col = 2)
hist(bif2[, 2], breaks = 20, xlab = "major-minor distance", main = "",xlim=c(0.002,0.01),ylim=c(0,5000))
abline(v = obs2[2], col = 2)

dev.off()
plot(density(mig2[, 2]))
lines(density(bif2[, 2]),col="red")
lines(density(hyb2[,2]),col="gray")
abline(v = obs2[2], col = 2)

models2<-rbind(bif2,mig2,hyb2)
models2<-data.frame(models2)
nrow(mig2)

index<-c(rep("bif",1000),
         rep("mig",1000),
         rep("hyb",1000))

x<-postpr(target=data.frame(obs2),index=index,sumstat=models2,tol=0.01, method="rejection")
summary(x)

x<-abc(target=data.frame(obs2),sumstat=hyb2,param=hyb[,1:7],tol=0.1,method="neuralnet")
plot(x,param=mig[,1:7])
summary(x)

x<-cv4abc(param=mig[,1:7],sumstat=mig2,nval = 100,tols=c(0.1),method="neuralnet")
plot(x)
