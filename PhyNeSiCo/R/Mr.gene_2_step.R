
sim.mr.gene<-function(tree, 
                      subst.rate, 
                      gamma.shape,
                      gamma.rate, 
                      number_loci, 
                      alleles, 
                      clock.like=T,
                      nsims){
table<-NULL
for(i in 1:nsims){
  result<-Mr.Gene2(tree, subst.rate=subst.rate, gamma.shape=gamma.shape, gamma.rate=gamma.rate, number_loci=number_loci, alleles=alleles, clock.like=clock.like)
  x<-cbind(result$sim_Principal_eigen, result$sim_asymmetry, result$sim_peakedness1, result$sim_peakedness2)
  y<-t(as.matrix(result$thetas))
  x<-cbind(y,x)
  table<-rbind(table,x)
  print(i)
}
colnames(table)<-c(paste("theta",c(1:length(result$thetas)),sep="_"),paste("eigen",c(1:ncol(result$sim_Principal_eigen)),sep="_"),"asymmetry","peakedness1","peakedness2")
return(table)
}

get.panda<-function(trees, method="standard"){
  results<-list()
  principal_eigenvalue<-NULL
  asymmetry<-NULL
  peakedness1<-NULL
  peakedness2<-NULL
  eigenvalues<-NULL
  for(i in 1:length(trees)){
    x<-spectR(trees[[i]], method=method)
    principal_eigenvalue<-rbind(principal_eigenvalue,x$principal_eigenvalue)
    asymmetry<-rbind(asymmetry,x$asymmetry)
    peakedness1<-rbind(peakedness1,x$peakedness1)
    peakedness2<-rbind(peakedness2,x$peakedness2)
    eigenvalues<-rbind(eigenvalues,x$eigenvalues)
   # print(i)
  }
  results$principal_eigenvalue<-principal_eigenvalue
  results$asymmetry<-asymmetry
  results$peakedness1<-peakedness1
  results$peakedness2<-peakedness2
  results$eigenvalues<-eigenvalues
  return(results)
}

table<-NULL
for(j in 1:100){

  spname <- species.name(write.tree(tree))
  seq <- rep(alleles, length(tree$tip.label))
  write.tree(tree)->tree2
  tree2<-put.thetas.on.nwk(tree2)
  mess.with.thetas(tree2, gamma.shape=2, gamma.rate=200)->tree_test
  out<-sim.coaltree.sp.mu(sptree=tree_test$tree, spname=spname, seq=seq, numgenetree=304,method="gamma",alpha=subst.rate)$gt
  sim.tre<-read.tree(text=out)
 
  xx<-NULL
  yy<-NULL
  for(i in 1:length(sim.tre)){
    x<-treedist(treeA,sim.tre[[i]])
    y<-treedist(treeB,sim.tre[[i]])
    xx<-rbind(xx,x)
    yy<-rbind(yy,y)
  }
  
  
  zz<-NULL
  for(ii in 1:nrow(xx)){
    z<-which(c(sum(xx[ii,c(1,3)]),sum(yy[ii,c(1,3)]))==min(sum(xx[ii,c(1,3)]),sum(yy[ii,c(1,3)])))
    if(length(z)>1){
      z<-0  
    }
    zz<-c(zz,z)
  }
  
  A<-length(grep(1,zz))/304
  B<-length(grep(2,zz))/304
  none<-length(grep(0,zz))/304
  
  
  
  A.2<-NULL
  B.2<-NULL
  for(i in 1:length(sim.tre)){
    a.2<-gsi(sim.tre[[i]],group1)
    b.2<-gsi(sim.tre[[i]],group2)
    A.2<-c(A.2,a.2)
    B.2<-c(B.2,b.2)
  }
  GSI<-cbind(A.2,B.2)
  y<-NULL
  for(i in 1:nrow(GSI)){
    x<-max(GSI[i,])
    y<-c(y,which(GSI[i,]==x))
  }
  GSI.a<-length(which(y==1))/nrow(GSI)
  GSI.b<-length(which(y==2))/nrow(GSI)
  
  table<-rbind(table,cbind(A,B,none,GSI.a,GSI.b))
  print(j)
  boxplot(table)
}

hist(table[,1])
abline(v=observed[1],col=2)
quantile(table[,1],probs=c(0.05,0.95))
observed[1]

hist(table[,2])
abline(v=observed[2],col=2)
quantile(table[,2],probs=c(0.05,0.95))
observed[2]

hist(table[,4])
abline(v=observed[4],col=2)
quantile(table[,4],probs=c(0.05,0.95))
observed[4]

hist(table[,5])
abline(v=observed[5],col=2)
quantile(table[,5],probs=c(0.05,0.95))
observed[5]

########################################################################
#####################################################################
#####################################################################3
index<-c(rep("tree1",nrow(table1)))
,rep("tree2",nrow(table2)))
table<-rbind(table1,table2)

table<-as.data.frame(table)

mod.prob<-abc(target=observed[1,1:10], param=table1[,1:45], sumstat=table1[,46:55],tol=0.1, sizenet=10,numnet=20,method="neuralnet", maxit=5000, trace=F)
plot(mod.prob, param=table1[,1:45])

summary(mod.prob)

fit<-NULL
for(i in 1:nrow(observed)){
x<-gfit(target=observed[i,c(1:10)],sumstat=table1[,c(46:55)],nb.replicate=100, trace=F,statistic = median)
plot(x)
x<-summary(x)
fit<-rbind(fit,x$pvalue)
print(c(i,x$pvalue))
}

fit2<-fit

length(which(fit<=0.05))/304

fit20summary(mod.prob)

mean(table1[5,1:45])

probs<-NULL
for(i in 1:nrow(observed)){
mod.prob<-abc(observed[1,1:10], index=index, table[,46:55],tol=0.01, sizenet=10,numnet=20,method="neuralnet", maxit=5000, trace=F)
x<-names(sort(mod.prob$pred, decreasing=T)[1])
x<-c(x,sort(mod.prob$pred, decreasing=T)[1])
probs<-rbind(probs,x)
print(i)
}

plot(density(as.numeric(probs[grep("tree2",probs[,1]),2])))
plot(density(as.numeric(probs[grep("tree1",probs[,1]),2])))

length(grep("tree2",probs))/nrow(probs)

x<-randomForest(theta_1 ~ ., data=table1[,c(1,46:ncol(table1))], na.action = na.omit, do.trace=T,ntree=2000, importance=T, proximity=T)
print(x)
round(importance(x), 2)
pred<-?predict(x,observed[1,])

colnames(table1[,c(1,46:ncol(table1))])
