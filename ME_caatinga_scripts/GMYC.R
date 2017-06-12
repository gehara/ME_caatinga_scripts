clade.cluster<-function(path.to.input,path.to.output,write.fastas,clust.prob,
                        min.samp,ntrees,mcmc,burnin){
  library(bGMYC)
  library(ape)
  #library(phytools)
  setwd(path.to.input)
  folders<-list.files()
  output<-list(NULL,NULL,NULL)
  names(output)<-c("result.single","result.multi","result.prob")
  result.singletree<-list()
  result.multitree<-list()
  result.prob<-list()
  for(i in 1:length(folders)){
    setwd(paste("./",folders[i],sep=""))
    files<-list.files()
    x<-grep("trees",files)
    trees<-read.nexus(files[x])
    trees<-trees[(length(trees)-8000):length(trees)]
    #### randomly resolve polytomies
   # for(i in 1:length(trees)){
    # tt<-is.binary.tree(trees[[i]])
     # if(tt==FALSE){
      #  trees[[i]] <- plot(multi2di(trees[[i]]))
      # }
    #}
    result<-bgmyc.singlephy(trees[[1]], mcmc=mcmc, burnin=1, thinning=10, t1=1, t2=10, start=c(1,1,2))
    result.singletree[[i]]<-result
    plot(result)
    pdf(file="result1.pdf",paper="a4r")
    plot(result)
    dev.off()
    trees<-sample(trees,ntrees)
    result.multi<-bgmyc.multiphylo(trees, mcmc=mcmc, burnin=burnin, thinning=100,start=c(1,1,15))
    pdf(file="result.multitree.pdf",paper="a4r")
    plot(result.multi)
    dev.off()
    plot(result.multi)
    result.multitree[[i]]<-result.multi
    result.spec<-bgmyc.spec(result.multi)
    result.probmat<-spec.probmat(result.multi)
    pdf(file="heatmap.pdf",paper="a4r")
    plot(result.probmat,trees[[1]])
    dev.off()
    plot(result.probmat,trees[[1]])
    result.prob[[i]]<-result.probmat
    bgmycrates<-checkrates(result.multi)
   
    clusters<-bgmyc.point(result.probmat,clust.prob)
    y<-grep("fas",files)
    fasta<-read.dna(files[y], format="fasta")
    fastas<-list()
    for(i in 1:length(clusters)){
    x<-c(which(rownames(fasta) %in% clusters[[i]]))
    fastas[[i]]<-fasta[c(x),]
  }
   if(write.fastas==T){ 
  for(i in 1:length(fastas)){
    setwd(path.to.output)
    write.dna(fastas[[i]],file=paste(strsplit(files[y],".fas")[[1]],"_",i,".fas",sep=""),format ="fasta")
  }}
  setwd(path.to.input)
  }
  
  setwd(path.to.output)
  fas.clust<-list.files()
  fas.clust<-fas.clust[grep(".fas",fas.clust)]
  for(i in 1:length(fas.clust)){
    x<-read.dna(fas.clust[i],format = "fasta")
    if(length(x[,1])<min.samp){
      file.remove(fas.clust[i])
    }
  }
  
  output$result.single<-result.singletree
  output$result.multi<-result.multitree
  output$result.prob<-result.prob
  return(output)
}

result<-clade.cluster(path.to.input="~/Dropbox/AMNH/Cophylogeography/GMYC",
                      path.to.output="~/Dropbox/AMNH/Cophylogeography/",
                      write.fastas=T,
                      clust.prob=0.5,
                      min.samp = 15,ntrees=100,mcmc=50000,burnin=40000)




