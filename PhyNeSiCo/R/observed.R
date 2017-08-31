
get.observed<-function(tree.b,
                       path,
                       distance,
                       spectral,
                       momentum){
  setwd(path)
  x<-list.files()
  x<-x[grep(".fas",x)]
  tab.tab<-NULL
for(i in 1:length(x)){
  
  fas<-read.dna(x[i],format="fasta")
  fas<-fas[match(tree.b[[1]],rownames(fas)),]
  length<-ncol(fas)
  fas<-fas2ms2(fas)
  fas<-ms.to.DNAbin(fas,length)
  
  if(distance==T){
    d<-dist.dna(fas, model="raw")
    sim.t<-as.vector(d)
  }

    if(spectral==T){
    d<-dist.dna(fas, model="raw")
    upgma.tree<-upgma(d)
    upgma.tree$edge.length<-upgma.tree$edge.length*1000
    eigen<-spectR(upgma.tree,"standard")
    sim.t<-c(eigen$principal_eigenvalue,eigen$asymmetry,eigen$eigenvalues)
    }
  
  tab.tab<-rbind(tab.tab,sim.t)
  print(i)
}
  if(momentum==T){
  observed<-c(apply(tab.tab,2,mean),apply(tab.tab,2,var),apply(tab.tab,2,kurtosis),apply(tab.tab,2,skewness))
  } else {observed<-tab.tab}
  
  return(observed)
}




