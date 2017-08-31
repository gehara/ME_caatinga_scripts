setwd("~/Dropbox/AMNH/Lampropeltis/Lampropeltis_TOL/Code")
source("network.sim.R")
astral.tree<-"(((((zonata:0.2694659048,multifasciata:0.2694659048):0.2699155957,(knoblochi:0.2694549107,pyromelana:0.2694549107):0.2699265898):0.269976768,alterna:0.8093582685):0.8102639922,((mexicana:0.4045359888,ruthveni:0.4045359888):0.4049105953,webbi:0.8094465841):0.8101756765):0.8098547822,((abnorma:0.6071502151,micropholis:0.6071502151):0.6073908336,polyzona:1.214541049):1.214935994);"
beast.tree<-"(((((zonata:1.534640091,multifasciata:1.534640091):1.423238548,(knoblochi:1.634034342,pyromelana:1.634034342):1.323844296):1.596522174,alterna:4.554400812):1.244838359,((mexicana:1.539058445,ruthveni:1.539058445):1.353206132,webbi:2.892264576):2.906974595):1.203233603,((abnorma:2.413413791,micropholis:2.413413791):2.046516169,polyzona:4.459929961):2.542542814);"
#par(mfrow=c(1,2))
#plot(read.tree(text=astral.tree))
#axisPhylo()
#plot(read.tree(text=beast.tree),direction = "leftwards")
#axisPhylo()

#get.tree.info(astral.tree)
#get.tree.info(beast.tree)
#2/2.4294770429
#5/2.4294770429
#2/7.002472775
#5/7.002472775

#setwd("/home/marcelo/Dropbox/AMNH/Lampropeltis/gene_fasta_11taxa/")
#observed.dist<-get.observed(tree.b=get.tree.info(astral.tree),path=getwd(),spectral=F, distance=T, momentum=T)
#observed.dist<-t(data.frame(observed.dist))
#observed.spec<-get.observed(tree.b=get.tree.info(astral.tree),path=getwd(),spectral=T, distance=F, momentum=F)

#hist(observed.dist[,1], breaks = 50)

#### Beast
### beast plot
#par(mfrow=c(1,2))
#plot(read.tree(text=get.tree.info(beast.tree)[[2]][1]))
#axisPhylo()
#plot(read.tree(text=beast.tree),direction = "leftwards")
#axisPhylo()

### beast.priors
#hib.priors<-c(2.9,4.4,4.4,0.5) # time of split: LOWER,UPPER; minor node UPPER; minor porportion;
#hib.clade<-c(6,7,8)
#major.sister<-c(1,2,3,4,5)
#minor.sister<-11

########## Astral
## astral plot
par(mfrow=c(1,2))
plot(read.tree(text=get.tree.info(astral.tree)[[2]][1]))
plot(read.tree(text=astral.tree),direction = "leftwards")
#axisPhylo()
### Astral priors
#"4 6"   "0.8094465841"
#"1 3"   "1.214541049" 
#"6 11"  "1.6196222606"
hib.priors<-c(0.81,1.2,1.2,0.28)
hib.clade<-c(6,7,8)
major.sister<-c(1,2,3,4,5)
minor.sister<-11

# segsites: if true calculate a tree with the simualted segregating sites.
# coaltrees: use the simulated coalescent tree
# distance: get the distances as output
# spectral: get spectral eigenvalues from Panda package
#getwd()

#write.table(bif.astral,"bif_astral_distance.txt",quote=F, row.names = F,col.names = F, sep="\t")
#write.table(bif.beast,"bif_beast_distance.txt",quote=F, row.names = F,col.names = F, sep="\t")
#write.table(net.beast,"net_beast_distance.txt",quote=F, row.names = F,col.names = F, sep="\t")
#write.table(net.astral,"net_astral_distance.txt",quote=F, row.names = F,col.names = F, sep="\t")

setwd("~/Dropbox/AMNH/Lampropeltis/Lampropeltis_TOL/Code")

for(i in 1:80){
  
  print(i)
  bif.beast<-sim.sp.tree(tree=beast.tree, Ne.prior = c(50000,1000000),time.modif=c(0.1,1),
                          nsims=2,bp=c(1489,285),mi=c(1e-10,1e-9),gen.time=2,
                          biffurcating = T,ntrees=304,segsites=T,coaltrees = F,
                          time.scalar=1000000, distance=T,spectral=F)
  
  print(i)
  x<-sim.sp.tree(tree=astral.tree, Ne.prior = c(50000,1000000),time.modif=c(0.8,2.1),
                          nsims=2,bp=c(1489,285),mi=c(1e-10,1e-9),gen.time=2,
                          biffurcating = T,ntrees=304,segsites=T,coaltrees = F,
                          time.scalar=1000000, distance=T,spectral=F)
  
  print(i)
  hib.priors<-c(2.9,4.4,4.4,0.5)
  net.beast<-sim.sp.tree(tree=beast.tree, Ne.prior = c(50000,1000000),time.modif=c(0.1,1),
                         hib.clade = hib.clade, hib.priors = hib.priors, major.sister = major.sister, 
                         minor.sister = minor.sister,
                         nsims=100,bp=c(1489,285),mi=c(1e-10,1e-9),gen.time=2,
                         biffurcating = F,ntrees=304,segsites=T,coaltrees = F,
                         time.scalar=1000000, distance=T,spectral=F)
  
  
 # mig.beast<-sim.sp.tree(tree=beast.tree, Ne.prior = c(500000,200000),time.modif=c(0.1,1),
#                         hib.clade = hib.clade, hib.priors = hib.priors, major.sister = major.sister, 
#                         minor.sister = minor.sister,
#                         nsims=100,bp=c(1489,285),mi=c(4e-10,2e-10),gen.time=2,
#                         biffurcating = F,ntrees=304,segsites=T,coaltrees = F,
#                         time.scalar=1000000, distance=T,spectral=F, migration = T, mig=c(0.1,1))
  print(i)
  hib.priors<-c(0.81,1.2,1.2,0.5)
  net.astral<-sim.sp.tree(tree=astral.tree, Ne.prior = c(50000,1000000),time.modif=c(0.8,2.1),
                         hib.clade = hib.clade, hib.priors = hib.priors, major.sister = major.sister, 
                        minor.sister = minor.sister,
                       nsims=100,bp=c(1489,285),mi=c(1e-10,1e-9),gen.time=2,
                      biffurcating = F,ntrees=304,segsites=T,coaltrees = F,
                     time.scalar=1000000, distance=T,spectral=F)
  
  write.table(bif.astral,"bif_astral_distance.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  write.table(bif.beast,"bif_beast_distance.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  write.table(net.beast,"net_beast_distance.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  write.table(net.astral,"net_astral_distance.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  
}
stop()

bif.beast<-read.table("bif_beast_distance.txt")
bif.astral<-read.table("bif_astral_distance.txt")
net.beast<-read.table("net_beast_distance.txt")
net.astral<-read.table("net_astral_distance.txt")




index<-c(rep("astral",nrow(bif.astral)),rep("beast",nrow(bif.beast)))
models<-rbind(bif.astral[,4:ncol(bif.astral)],bif.beast[,2:ncol(bif.beast)])

index<-c(rep("bif.beast",nrow(bif.beast)),rep("net.beast",nrow(net.beast)),
         rep("bif.astral",nrow(bif.astral)),rep("net.astral",nrow(net.astral)))

models<-rbind(bif.beast[,6:ncol(bif.astral)],net.beast[,6:ncol(bif.astral)],
              bif.astral[,6:ncol(bif.astral)],net.astral[,6:ncol(bif.astral)])

models<-rbind(net.astral[,6:ncol(bif.astral)])

ncol(models)
ncol(observed.dist)
head(bif.astral)
hist(models[2,])


library(ggplot2)
colnames(observed.dist)<-NULL
colnames(models)->colnames(observed.dist)

xx<-rbind(models[,1:55],observed.dist[1:55])
xx<-rbind(models[,],observed.dist[])
data<-c(index,rep("observed",nrow(observed.dist)))

PCA<-prcomp(xx, center = T, scale. = T, retx=T)
plot(PCA)
scores <- data.frame(PCA$x[,1:9])
ggplot(scores, aes(x=PC1, y=PC2))+
  theme(legend.text = element_text(size = 16, face = "bold"))+
  geom_point(aes(colour=data, size=data, shape=data))+
  scale_size_manual(values=c(2,2,2,2,10))+
  scale_colour_brewer(palette="Set1")
dev.off()


w<-cv4postpr(index,sumstat=scores[1:(nrow(scores)-nrow(observed.dist)),],tols=0.01, 
             method="rejection", nval=100, maxit = 2000)
summary(w)
plot(w)

w<-cv4postpr(index,sumstat=models,tols=0.0005*150000, 
             method="rejection", nval=50)
summary(w)
plot(w)

ncol(models)
ncol(observed.dist[1,])
nrow(scores)
probs<-NULL

  x<-postpr(target=scores[nrow(scores),],index=data[1:(nrow(scores)-nrow(observed.dist))],sumstat=scores[1:(nrow(scores)-nrow(observed.dist)),],tol=0.01, method="neuralnet")
  y<-abc(scores[nrow(scores),],net.astral[,1:5],scores[1:(nrow(scores)-nrow(observed.dist)),],tol=0.02, method="neuralnet")
  summary(y)
  plot(y, param=net.astral[,1:5])


  goodness<-gfit(observed.dist,net.astral[,6:ncol(net.astral)],tol=0.05,nb.replicate = 100)
  summary(goodness)
  plot(goodness)
 x<-cv4abc(net.astral[,1:5],scores[1:(nrow(scores)-nrow(observed.dist)),], nval=100,
           method="neuralnet",tols=c(0.1,0.01))  
plot(x)  
  
  
xx<-NULL
for(i in 1:nrow(probs)){
  x<-which(probs[i,]==max(probs[i,]))[1]
  xx<-c(xx,x)
}

length(which(xx==1))/length(xx)
length(which(xx==2))/length(xx)


########################
################################
### beast.priors
hib.priors<-c(2.9,4.4,4.4,0.5) # time of split: LOWER,UPPER; minor node UPPER; minor porportion;
hib.clade<-c(9,10,11)
major.sister<-c(4,5,6,7,8)
minor.sister<-3

getwd()
setwd("~/Dropbox/AMNH/Lampropeltis/Lampropeltis_TOL/Code")
write.table(net,"net_beast_spectral.txt",quote=F, row.names = F,col.names = F, sep="\t")

for(i in 1:32){
  
  #bif<-sim.sp.tree(tree=beast.tree, Ne.prior = c(100000,1000000),time.modif=c(0.3,1),
  #                       nsims=1000,bp=c(1489,285),mi=c(4e-10,1e-10),gen.time=2,
  #                       biffurcating = T,ntrees=1,segsites=T,coaltrees = F,
  #                       time.scalar=1000000, distance=T,spectral=F)
  
  net<-sim.sp.tree(tree=beast.tree, Ne.prior = c(100000,1000000),time.modif=c(0.3,1),
                   hib.clade = hib.clade, hib.priors = hib.priors, major.sister = major.sister, 
                   minor.sister = minor.sister,
                   nsims=1000,bp=c(1489,285),mi=c(4e-10,1e-10),gen.time=2,
                   biffurcating = F,ntrees=1,segsites=T,coaltrees = F,
                   time.scalar=1000000, distance=F,spectral=T)
  
  
  write.table(net,"net_beast_spectral.txt",quote=F, row.names = F,col.names = F, sep="\t", append = T)
  #write.table(bif,"bif_beast_distance.txt",quote=F, row.names = F,col.names = F, sep="\t",append = T)
  
}

nrow(bif)
nrow(net)
bif<-read.table("bif_beast_spectral.txt")
net<-read.table("net_beast_spectral.txt")

index<-c(rep("bifurcate",nrow(bif)),rep("network",nrow(net)))
models<-rbind(bif[,2:ncol(bif)],net[,2:ncol(net)])

plot(bif[,2]~bif[,3])
points(net[,2]~net[,3], col="red")

########### PCA
###########
#colnames(observed)<-colnames(models)
library(ggplot2)
xx<-rbind(models,observed.spec)
data<-c(index,rep("observed",nrow(observed.spec)))

PCA<-prcomp(xx, center = T, scale. = T, retx=T)
plot(PCA)
scores <- data.frame(PCA$x[,1:10])
ggplot(scores, aes(x=PC1, y=PC2))+
  theme(legend.text = element_text(size = 16, face = "bold"))+
  geom_point(aes(colour=data, size=data, shape=data))+
  scale_size_manual(values=c(2,2,1))+
  scale_colour_brewer(palette="Set1")
dev.off()
############

##########
goodness<-gfit(observed[1,],net[,2:ncol(net)],tol=0.01,nb.replicate = 100)
summary(goodness)
plot(goodness)

goodness<-gfit(observed.dist[2,],bif.astral[,2:ncol(bif.astral)],tol=0.01,nb.replicate = 50)
summary(goodness)
plot(goodness)

goodness<-gfit(observed.dist[2,],bif.beast[,2:ncol(bif.beast)],tol=0.01,nb.replicate = 50)
summary(goodness)
plot(goodness)

w<-cv4postpr(index,sumstat=scores[1:(nrow(scores)-nrow(observed.dist)),],tols=0.001, 
             method="rejection", nval=50)
summary(w)
plot(w)

w<-cv4postpr(index,sumstat=models,tols=0.001, 
             method="rejection", nval=100)
summary(w)
plot(w)
nrow(scores)
probs<-NULL
for(i in 1:nrow(observed.spec)){
  
  x<-postpr(scores[(64000+i),],index,scores[1:64000,],tol=0.001, method="rejection")
  
  x<-summary(x)
  p<-x$Prob
  probs<-rbind(probs,p)
  
}
xx<-NULL
for(i in 1:nrow(probs)){
  x<-which(probs[i,]==max(probs[i,]))[1]
  xx<-c(xx,x)
}

length(which(xx==1))/length(xx)
length(which(xx==2))/length(xx)


length(grep(3,xx))/304

x<-postpr(observed,index,models,tol=0.1, method="rejection")
summary(x)

post<-abc(target=observed[],sumstat = bif[,2:ncol(bif)], param=bif[,1], tol=0.01,method="neuralnet")
plot(density(post$adj.values))
summary(post)
plot(post, param=bif[,1])

#######################################
######################################
#####################################33
library(help=neuralnet)

model<-c(rep("0",nrow(bif)),rep("1",nrow(net)))
models<-rbind(bif,net)

models<-data.frame(models)
head(models)
maxs <- apply(models, 2, max) 
mins <- apply(models, 2, min)
scaled <- as.data.frame(scale(models, center = mins, scale = maxs - mins))

scaled<-cbind(model,scaled)
head(scaled)

index <- sample(1:nrow(scaled),round(0.75*nrow(scaled)))

train_ <- scaled[index,]
test_ <- scaled[-index,]


library(help=neuralnet)
n <- names(train_)
f <- as.formula(paste("model ~", paste(n[!n %in% "model"], collapse = " + ")))
nn <- neuralnet(f,data=train_,lifesign="full", hidden=15, err.fct="sse", 
                linear.output=FALSE)
summary(nn)
plot(nn)
prediction(nn)

predict(nn, data=Test,type = "class")

nn <- neuralnet(f,data=train_,hidden=c(5,3),linear.output=T)




### calculate observed
setwd("~/Dropbox/AMNH/Lampropeltis/TICR/my-genes-mb/trees.nwk")
x<-list.files()
obs.tre<-read.tree(x[307])
group<-c("alterna","zonata","multifasciata","pyromelana","knoblochi","webbi",
         "ruthveni","mexicana","polyzona","abnorma","micropholis")
out<-obs.tre[[1]]$tip.label[-match(group,obs.tre[[1]]$tip.label)]
for(i in 1:length(obs.tre)){
  obs.tre[[i]]<-drop.tip(obs.tre[[i]],out)
}

treeA
observed<-NULL
for(i in 1:length(obs.tre)){
  x<-treedist(treeA,obs.tre[[i]])
  observed<-rbind(observed,x)
}

###############

boxplot(simulated[,1],observed[,1])
obs<-c(mean(observed[,1]),var(observed[,1]),kurtosis(observed[,1]),skewness(observed[,1]))
obs<-c(obs,mean(observed[,3]),var(observed[,3]),kurtosis(observed[,3]),skewness(observed[,3]))

plot(density(simulated[,1]))
lines(density(observed[,1]),col="red")


plot.sim.obs(bif,obs)
boxplot(cbind(simulated[,4],simulated.net[,4]))




