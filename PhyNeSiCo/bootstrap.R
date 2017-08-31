#combine all bootstrap replicates from diferent folders in a single list of trees
comb.replicates<-function(path, write.file=T){
  setwd(path)
  # search for outputs
  y<-NULL
  for(i in 1:100){
    setwd(paste("./",i,sep=""))
    x<-list.files(pattern="bestNet.out")
    if(length(x)==1){
      y<-c(y,i)
    }
    setwd("../")
  }
  i<-NULL
  trees<-NULL
  
  # combine outputs in a sigle object
  wd<-getwd()
  for(i in y){
    setwd(paste("./",i,sep=""))
    x<-readLines(con="bestNet.out",n=3)[3]
    trees<-c(trees,x)
    setwd(wd)
  }
  trees<-gsub("Dendroscope:","",trees)
  trees<-gsub(" ","",trees)
  writeLines(trees,con="replictrees.txt")
  return(y)
}
# exclude characters and branch lengths
getrid<-function(x){
  z<-c(0:9,".",":")
  for(i in z){
    x<-gsub(i,"",x,fixed=T)
  }
  return(x)
}
## change taxon names to numbers
taxons.to.numbers<-function(tree,namelist=F){
  
  if(namelist[1]==F){
    x<-tree  
  x<-gsub(")","",x,fixed=T)
  x<-gsub("(","",x,fixed=T)
  x<-gsub("#H","",x,fixed=T)
  x<-gsub(";","",x,fixed=T)
  x<-strsplit(x,",",fixed=T)[[1]]
  x<-x[x!=""]
  } else { x<-namelist}
  
  for(i in 1:length(x)){
    tree<-gsub(x[i],i,tree)
  }
  output<-list(x,tree)
  return(output)
}
### get the taxa in the tree
get.taxa<-function(x){
  x<-gsub("(","",x,fixed=T)
  x<-gsub(")","",x,fixed=T)
  x<-gsub(";","",x,fixed=T)
  x<-gsub("#H","",x,fixed=T)
  x<-strsplit(x,",",fixed=T,useBytes=T)[[1]]
  x<-sort(as.numeric(x))
  return(x)
}
#### get the hybrid node from tree
get.hybrid<-function(x){
  where<-1+regexpr(")#H",x)[1]
  x<-strsplit(x,"",fixed=T)[[1]]
  x<-x[-c(where:length(x))]
  zz<-0
  for(i in length(x):1){
    
    if(x[i]==")"){
      zz<-zz+1
    }
    if (x[i]=="("){
      zz<-zz-1
    }
    if(zz==0)break
    i<-i-1
  }
  hybrid.tree<-paste(x[i:length(x)],collapse="")
  hybrids<-get.taxa(hybrid.tree)
  return(list(hybrid.tree,hybrids))
}
## get parents from tree
get.parents<-function(bestNet,hybrid){
  x<-bestNet
  x<-gsub(hybrid,"",x,fixed=T)
  
  where<-4+regexpr(",#H",x)[1]
  x<-strsplit(x,"",fixed=T)[[1]]
  x<-x[-c(where:length(x))]
  zz<-0
  for(i in length(x):1){
    
    if(x[i]==")"){
      zz<-zz+1
    }
    if (x[i]=="("){
      zz<-zz-1
    }
    if(zz==0)break
    i<-i-1
  }
  first.parent<-paste(x[i:length(x)],collapse="")
  
  x<-bestNet
  x<-gsub(hybrid,"",x,fixed=T)
  x<-gsub(first.parent,"",x, fixed=T)
  
  where<-4+regexpr(",#H",x,fixed=T)[1]
  x<-strsplit(x,"",fixed=T)[[1]]
  x<-x[-c(where:length(x))]
  zz<-0
  for(i in length(x):1){
    if(x[i]==")"){
      zz<-zz+1
    }
    if (x[i]=="("){
      zz<-zz-1
    }
    if(zz==0)break
    i<-i-1
  }
  secound.parent<-paste(x[i:length(x)],collapse="")
  #
  firstparent<-cbind(first.parent,paste(get.taxa(first.parent),collapse=" "))
  secoundparent<-cbind(secound.parent,paste(get.taxa(secound.parent),collapse=" "))
  #
  return(rbind(firstparent,secoundparent))
}
### get hybrid node support
get.hybrid.support<-function(trees,BNet.hybrids){
  BoRep.hybrids<-NULL
  for(i in 1:length(trees)){
    BoRep.hybrids[i]<-paste(get.hybrid(trees[i])[[2]],collapse=" ")
  }
  hybrid.support<-length(which(BoRep.hybrids==paste(BNet.hybrids[[2]],collapse=" ")))/length(trees)
  hybrid.support<-round(hybrid.support,2)
  return(hybrid.support)
}
## get parents support
get.parents.support<-function(trees,BNet.parents){
  BoRep.parents<-matrix(ncol=2,nrow =length(trees))
  for(i in 1:length(trees)){
    hybrid<-paste(get.hybrid(trees[i])[[1]],collapse=" ")
    BoRep.parents[i,]<-get.parents(bestNet=trees[i],hybrid)[,2]
  }
  parent1<-length(which(BoRep.parents==BNet.parents[1,2]))/length(trees)
  parent2<-length(which(BoRep.parents==BNet.parents[2,2]))/length(trees)
  parent1<-round(parent1,2)
  parent2<-round(parent2,2)
  BNet.parents<-cbind(BNet.parents,c(parent1,parent2))
  colnames(BNet.parents)<-c("topology","parents","support")
  return(BNet.parents)
}
## get reticulation support
get.reticulation.support<-function(trees,BNet.parents,BNet.hybrids){
  
  BoRep<-matrix(ncol=3,nrow=length(trees))
  
  for(i in 1:length(trees)){
    hyb<-paste(get.hybrid(trees[i])[[1]],collapse=" ")
    BoRep[i,1]<-paste(get.hybrid(trees[i])[[2]],collapse=" ")
    BoRep[i,2:3]<-get.parents(bestNet=trees[i],hyb)[,2]
  }
  
  freq<-0
  for(i in 1:nrow(BoRep)){
    x<-length(which(BoRep[i,1]==paste(BNet.hybrids[[2]],collapse=" ")))
    y<-length(which(BoRep[i,2:3]==paste(BNet.parents[1,2],collapse=" ")))
    z<-length(which(BoRep[i,2:3]==paste(BNet.parents[2,2],collapse=" ")))
    if((x+z+y)==3){
      freq<-freq+1  
    }
  }
  ret.support<-freq/length(trees)
  ret.support<-round(ret.support,2)
  return(ret.support)
}

trees<-comb.replicates(path="~/Desktop/bootstrap/")
comb.replicates(path="~/Desktop/bootstrap2/")

trees<-readLines("~/Desktop/bootstrap/replictrees.txt")
trees2<-readLines("~/Desktop/bootstrap2/replictrees.txt")
trees<-c(trees,trees2)

bestNet<-readLines("~/Desktop/bestnet.txt")

bestNet<-getrid(bestNet)

bestNet<-taxons.to.numbers(bestNet)

trees<-getrid(trees)

for(i in 1:length(trees)){
trees[i]<-taxons.to.numbers(trees[i],namelist=bestNet[[1]])[[2]]
}
trees <- trees[trees!="NA"]

BNet.hybrids <- get.hybrid(bestNet[[2]])

Hyb.support <- get.hybrid.support(trees,BNet.hybrids)

BNet.parents <- get.parents(bestNet[[2]],BNet.hybrids[[1]])

BNet.parents <- get.parents.support(trees,BNet.parents)

ret.support <- get.reticulation.support(trees,BNet.parents,BNet.hybrids)

BoRep.parents<-matrix(ncol=3,nrow=length(trees)*2)
j<-1
for(i in 1:length(trees)){
  hybrid<-paste(get.hybrid(trees[i])[[1]],collapse=" ")
  x<-get.parents(bestNet=trees[i],hybrid)
  x<-get.parents.support(trees,x)
  BoRep.parents[j:(j+1),]<-x
  j<-j+2
}

BoRep.hyb<-matrix(ncol=3,nrow=length(trees))
for(i in 1:length(trees)){
  hybrid<-get.hybrid(trees[i])
  x<-get.hybrid.support(trees,hybrid)
  BoRep.hyb[i,1]<-hybrid[[1]]
  BoRep.hyb[i,2]<-paste(hybrid[[2]],collapse=" ")
  BoRep.hyb[i,3]<-x
}


BoRep.reticulation<-matrix(nc=4,nrow=length(trees))
for(i in 1:length(trees)){
  hybrid<-get.hybrid(trees[i])
  parents<-get.parents(bestNet=trees[i],hybrid)
  parents<-parents[order(parents[,2]),]
  x<-get.reticulation.support(trees,parents,hybrid)
  BoRep.reticulation[i,1]<-paste(hybrid[[2]], collapse=" ")
  BoRep.reticulation[i,2]<-paste(parents[1,2], collapse=" ")
  BoRep.reticulation[i,3]<-paste(parents[2,2], collapse=" ")
  BoRep.reticulation[i,4]<-x
}

x<-unique(BoRep.reticulation)
x<-x[rev(order(as.numeric(x[,4]))),]
colnames(x)<-c("hybrid","parent1","parent2","support")
write.table(x,file="reticulation.txt",quote=F,col.names = T,row.names = F, sep="\t")


leg<-data.frame(matrix(nc=2,nr=length(bestNet[[1]])))
leg[,1]<-c(1:length(bestNet[[1]]))
leg[,2]<-bestNet[[1]]
colnames(leg)<-c("Number","Species")

all.sup<-list("NULL","NULL","NULL")
names(all.sup)<-c("legend","Hybrids","parents")
x<-BoRep.parents[,2:3]
x<-unique(x)
x<-x[rev(order(x[,2])),]
colnames(x)<-c("species","support")
x<-data.frame(x)
all.sup$parents<-x

x<-unique(BoRep.hyb[,2:3])
x<-x[rev(order(x[,2])),]
colnames(x)<-c("species","support")
x<-data.frame(x)
all.sup$Hybrids<-x

all.sup$legend<-leg

write.table(all.sup$legend,file="legend.txt",quote=F,col.names = T,row.names = F, sep="\t")
write.table(all.sup$Hybrids,file="hybrids.txt",quote=F,col.names = T,row.names = F, sep = "\t")
write.table(all.sup$parents,file="parents.txt",quote=F,col.names = T,row.names = F, sep = "\t")
getwd()


