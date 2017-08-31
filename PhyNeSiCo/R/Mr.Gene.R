###simulate gene trees from species trees given edge.lengths (substitution rates) and theta and generate summary stats 
###tree=phylo. object, subst.rate=distribution of edge.lengths representing substitution lengths,gamma.shape=for theta generation,gamma.rate=for theta generation, number_loci=number of gene trees to simulate, alleles=alleles per taxon, assumed 1)
###By Burbrink and Gehara, uses the phybase engine. 


library(phybase)
library(RPANDA)
library(phangorn)



Mr.Gene2<-function(tree, subst.rate, gamma.shape,gamma.rate, number_loci, alleles, clock.like=T)
  
{
  put.thetas.on.nwk<-function(tree){
    x<-strsplit(tree,":")[[1]]
    for(i in 2:length(x)){
      y<-strsplit(x[i],"")[[1]]
      for(j in 1:length(y)){
        if(y[j] %in% c(",",")")){
          x[i]<-paste(paste(y[1:(j-1)],collapse=""),"#0.01",paste(y[j:length(y)],collapse=""),sep="")
          next
        } 
      }
    }
    x[length(x)]<-gsub(";","#0.01;",x[length(x)])
    x<-paste(x,collapse=":")
    return(x)
  }
  
  
  ### this is part of the main function bellow and
  # is used to get the thetas from the tree string 
  # ex of string input: "(((S1:0.005#0.01,S2:0.005#0.01):0.00025#.01,S3:0.00525#0.01):0.00025#0.01,S4:0.0055#0.01)#.01;"
  
  get.thetas<-function(tree){
    tree<-strsplit(tree,"#")[[1]]
    thetas<-NULL
    for(i in 2:length(tree)){
      
      x<-strsplit(tree[i],"")[[1]]
      
      for(j in 1:length(x)){
        
        if(x[j] %in% c(",",")",";")){
          t<-paste(x[1:(j-1)],collapse="")
        }
      }
      thetas<-c(thetas,t)
    }
    return(thetas)
  }
  
  # this function replace thetas with new sampled values from a gamma distribution
  # It uses the function above to get the thetas from the tree string.
  # You need a tree string as input and two gamma parameters (gamma.shape and gamma.rate)
  
  # It returns a list with the new tree string as the first element, and the new thetas as a second element.
  
  mess.with.thetas<-function(tree,gamma.shape,gamma.rate){
    
    thetas<-get.thetas(tree)
    new.thetas<-round(rgamma(length(thetas),gamma.shape,gamma.rate),4)
    tree<-strsplit(tree,"#")[[1]]
    for(i in 1:length(thetas)){
      tree[i+1]<-gsub(thetas[i],new.thetas[i],tree[i+1])
    }
    tree<-paste(tree,collapse="#")
    output<-list()
    output$tree<-tree
    output$thetas<-new.thetas
    return(output)
  }
  
  
  
  #species name
  spname <- species.name(write.tree(tree))
  seq <- rep(alleles, length(tree$tip.label))
  
  #print("simulating trees")
  results<-list()
  dd<-NULL
  #for(i in 1:1000){
  if(clock.like==F){
    write.tree(tree)->tree2
    tree2<-put.thetas.on.nwk(tree2)
    mess.with.thetas(tree2, gamma.shape, gamma.rate)->tree_test
    results$species.tree<-tree_test$tree
    results$thetas<-tree_test$thetas
    
  out<-sim.coaltree.sp.mu(sptree=tree_test$tree, spname=spname, seq=seq, numgenetree=304,method="gamma",alpha=subst.rate)$gt
  }
  
  if(clock.like==T){
    ###get nodematrix
    nodematrix <- read.tree.nodes(write.tree(tree), spname)$nodes
    nodematrix[,5]<-rgamma(nrow(nodematrix),gamma.shape,gamma.rate)
    ##rootnode
    rootnode <- rootoftree(nodematrix)
    
  out<-sim.coaltree.sp(rootnode=rootoftree(nodematrix),nodematrix=nodematrix,nspecies=length(tree$tip.label),seq=seq,name=spname)$gt
  }
  #d<-read.tree(text=paste(out))
  dd<-c(dd,out)
  #print(i)
  #}
  
  simulated<-lapply(c(1:length(dd)), function(x) read.tree(file="",text=dd[x]))
  class(simulated)<-"multiPhylo"
  results$gene_trees<-simulated
  
  
  #get spectral densities
  #print("getting spectral densities")
  
  sim_P_eigen<-NULL
  sim_P_asymmetry<-NULL
  sim_P_peakedness1<-NULL
  sim_P_peakedness2<-NULL
  
  for(i in 1:length(simulated)){
    x<-spectR(simulated[[i]], "normal")
    sim_P_eigen<-rbind(sim_P_eigen,x$eigenvalue)
    sim_P_asymmetry<-c(sim_P_asymmetry,x$asymmetry)
    sim_P_peakedness1<-c(sim_P_peakedness1,x$peakedness1)
    sim_P_peakedness2<-c(sim_P_peakedness1,x$peakedness2)[2]
    #print(i)
  }
  
  
  results$sim_Principal_eigen<-sim_P_eigen
  results$sim_asymmetry<-sim_P_asymmetry
  results$sim_peakedness1<-sim_P_peakedness1
  results$sim_peakedness2<-sim_P_peakedness2
  
  #results$rf<-unlist(lapply(simulated, function(x) RF.dist(x, tree)))
  
  #print("done!")
  return(results)
 
}

randomT<-rcoal(length(spname), rooted = TRUE, tip.label = spname)
spname <- species.name(write.tree(tree1))
plot(randomT)
t<-write.tree(randomT)
randomT<-read.tree(text=t)

treeA<-"(Cemophora_coccinea:11.678229466235559,(((Lampropeltis_elapsoides:6.240913466415052,(Lampropeltis_annulata:2.1388501466619054,Lampropeltis_gentilis:2.1388501466619054):4.102063319753147):2.3046166642738006,(((Lampropeltis_abnorma:2.413413791353035,Lampropeltis_micropholis:2.413413791353035):2.046516169196659,Lampropeltis_polyzona:4.459929960549694):2.5425428142366115,((Lampropeltis_alterna:4.554400812401852,((Lampropeltis_zonata:1.5346400910659765,Lampropeltis_multifasciata:1.5346400910659765):1.4232385476943845,(Lampropeltis_knoblochi:1.6340343424803274,Lampropeltis_pyromelana:1.6340343424803274):1.3238442962800336):1.5965221736414925):1.2448383594795835,(Lampropeltis_webbi:2.892264576497898,(Lampropeltis_mexicana:1.539058444509318,Lampropeltis_ruthveni:1.539058444509318):1.3532061319885798):2.9069745953835375):1.20323360290487):1.5430573559025476):1.9194000035200798,(((Lampropeltis_nigra:2.623141453795378,Lampropeltis_getula:2.623141453795378):1.7038376207295447,(Lampropeltis_holbrooki:2.969604917867215,(Lampropeltis_californiae:1.6295385079701354,Lampropeltis_splendida:1.6295385079701354):1.3400664098970814):1.3573741566577078):2.3942359444784973,(Lampropeltis_extenuata:5.480392736726901,(Lampropeltis_rhombomaculata:2.3301760341467883,Lampropeltis_calligaster:2.3301760341467883):3.1502167025801127):1.2408222822765191):3.7437151152055126):1.2132993320266259);"
treeA<-read.tree(text=treeA)
treeB<-"(Cemophora_coccinea:11.646531622834841,(((Lampropeltis_elapsoides:6.4724541318431905,(Lampropeltis_annulata:2.223557678670028,Lampropeltis_gentilis:2.223557678670028):4.248896453173162):1.7351824579374346,(((Lampropeltis_abnorma:2.299446356630222,Lampropeltis_micropholis:2.299446356630222):3.015957032048095,(Lampropeltis_polyzona:4.494348796582882,(Lampropeltis_webbi:2.437111922702547,(Lampropeltis_mexicana:1.2095584711996956,Lampropeltis_ruthveni:1.2095584711996956):1.2275534515028532):2.057236873880335):0.8210545920954351):1.5441590364726077,(Lampropeltis_alterna:4.98585614680508,((Lampropeltis_zonata:1.6551188768199818,Lampropeltis_multifasciata:1.6551188768199818):1.498265321401723,(Lampropeltis_knoblochi:1.7104033936538254,Lampropeltis_pyromelana:1.7104033936538254):1.4429808045678794):1.8324719485833754):1.8737062783458445):1.3480741646297005):2.0683784200558453,(((Lampropeltis_nigra:2.7940973017001376,Lampropeltis_getula:2.7940973017001376):1.6979135085595853,(Lampropeltis_holbrooki:3.1538217772605943,(Lampropeltis_californiae:1.7761760182728263,Lampropeltis_splendida:1.7761760182728263):1.3776457589877698):1.3381890329991286):2.4331913359368373,(Lampropeltis_extenuata:5.663054476446628,(Lampropeltis_rhombomaculata:2.453455947146429,Lampropeltis_calligaster:2.453455947146429):3.209598529300199):1.2621476697499325):3.35081286363991):1.3705166129983706);"  
treeB<-read.tree(text=treeB)
plot(treeA)
plot(treeB)

treeA$edge.length<-treeA$edge.length/1000
#treeB$edge.length<-treeB$edge.length/1000

for(i in 1:length(treeA$tip.label)){
  treeA$tip.label[i]<-strsplit(treeA$tip.label,"_")[[i]][2]
}

for(i in 1:length(treeB$tip.label)){
  treeB$tip.label[i]<-strsplit(treeB$tip.label,"_")[[i]][2]
}

xx<-NULL
yy<-NULL
for(i in 1:length(tre)){
x<-treedist(treeA,tre[[i]])
y<-treedist(treeB,tre[[i]])
xx<-rbind(xx,x)
yy<-rbind(yy,y)
}

zz<-NULL
for(i in 1:nrow(xx)){
z<-which(c(sum(xx[i,c(1,3)]),sum(yy[i,c(1,3)]))==min(sum(xx[i,c(1,3)]),sum(yy[i,c(1,3)])))
if(length(z)>1){
z<-0  
}
zz<-c(zz,z)
}

a<-length(grep(1,zz))/304
b<-length(grep(2,zz))/304
n<-length(grep(0,zz))/304

observed<-c(a,b,n)

boxplot(xx[1,1],yy[1,1])

tre[[1]]$edge.length

plot(density(rgamma(1000,2,200)))

median(rgamma(1000,2,100))

table1<-sim.mr.gene(tree1, subst.rate=rnorm(1,1,0.1), 
                    gamma.shape=10, gamma.rate=1000, 
                    number_loci=1, alleles=1, clock.like=F,
                    nsims=1000)

table2<-sim.mr.gene(randomT, subst.rate=rnorm(1,1,0.1), 
                    gamma.shape=0.001, gamma.rate=1, 
                    number_loci=1, alleles=1, clock.like=F,
                    nsims=1000)

setwd("~/Dropbox/AMNH/Lampropeltis/TICR/my-genes-mb/trees.nwk")
library(ape)
x<-list.files()
tre<-read.tree(x[307])
treedist(tree1,tre[[1]])

observed<-get.panda(tre)
observed<-cbind(observed$sim_Principal_eigen,observed$sim_asymmetry, observed$sim_peakedness1, observed$sim_peakedness2)

colnames(observed)<-c(paste("eigen",c(1:ncol(result$sim_Principal_eigen)),sep="_"),"asymmetry","peakedness1","peakedness2")



plot(table1[,46]~table1[,47])
points(observed[,1]~observed[,2], col="red")

points(table2[,46]~table2[,47], col="blue")

plot(observed[,1], col="red")
points(table2[,1])
