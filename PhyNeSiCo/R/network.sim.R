library(RPANDA)
library(phybase)
library(phyclust)
library(phangorn)
library(PipeMaster)
library(ape)
library(e1071)
library(abc)
library(msm)

##### get observed 
fas2ms2<-function (fas) 
{
  
  fas <- as.character(fas)
  
  bin<-NULL
  pos<-NULL
  gaps<-0
  for(i in 1:ncol(fas)){
    a<-length(grep("a",fas[,i]))
    c<-length(grep("c",fas[,i]))
    g<-length(grep("g",fas[,i]))
    t<-length(grep("t",fas[,i]))
    
    r<-length(grep("r",fas[,i]))
    y<-length(grep("y",fas[,i]))
    m<-length(grep("m",fas[,i]))
    k<-length(grep("k",fas[,i]))
    s<-length(grep("s",fas[,i]))
    w<-length(grep("w",fas[,i]))
    
    h<-length(grep("h",fas[,i]))
    b<-length(grep("b",fas[,i]))
    v<-length(grep("v",fas[,i]))
    d<-length(grep("d",fas[,i]))
    n<-length(grep("n",fas[,i]))
    gap<-length(grep("-",fas[,i]))
    if(gap>0){
      gaps<-gaps+1
      next
    } else if (n>0){next}
    if(!nrow(fas) %in% c(a,c,g,t,r,y,m,k,s,w,h,b,v,d)){
      bin<-cbind(bin,fas[,i])
      pos<-c(pos,i)
    }
  }
  pos<-pos/ncol(fas)
  if(!(is.null(bin))){
    for(j in 1:ncol(bin)){
      g<-length(grep("g",bin[,j]))/nrow(bin)
      a<-length(grep("a",bin[,j]))/nrow(bin)
      t<-length(grep("t",bin[,j]))/nrow(bin)
      c<-length(grep("c",bin[,j]))/nrow(bin)
      
      r<-length(grep("r",bin[,j]))/nrow(bin)
      y<-length(grep("y",bin[,j]))/nrow(bin)
      m<-length(grep("m",bin[,j]))/nrow(bin)
      k<-length(grep("k",bin[,j]))/nrow(bin)
      s<-length(grep("s",bin[,j]))/nrow(bin)
      w<-length(grep("w",bin[,j]))/nrow(bin)
      
      h<-length(grep("h",bin[,j]))/nrow(bin)
      b<-length(grep("b",bin[,j]))/nrow(bin)
      v<-length(grep("v",bin[,j]))/nrow(bin)
      d<-length(grep("d",bin[,j]))/nrow(bin)
      bases<-c(a,c,g,t,r,y,m,s,k,w,h,b,v,d)
      names(bases)<-c("a","c","g","t","r","y","m","s","k","w","h","b","v","d")
      base<-which.max(bases[1:4])
      bin[,j]<-gsub(names(base),0,bin[,j])
      bases<-bases[-base]
      for(i in 1:length(bases)){
        bin[,j]<-gsub(names(bases[i]),1,bin[,j])
      }
    }
    seqs<-NULL
    for(i in 1:nrow(bin)){
      seqs<-c(seqs,paste(bin[i,],collapse=""))
    }
    ss<-ncol(bin)
  }else{seqs<-NULL
  ss<-0}
  
  
    ms.out <- paste("ms", nrow(fas), 1)
    ms.out <- c(ms.out, paste("//"))
    ms.out <- c(ms.out, paste("segsites:", 
                              ncol(bin)))
    ms.out <- c(ms.out, paste("positions:   ", 
                              paste(pos, collapse = "    ")))
    ms.out <- c(ms.out, paste(seqs, sep = "\n"))

   return(list(ms.out,gaps))

}


get.observed<-function(tree.b,
                       path,
                       distance,
                       spectral,
                       moments){
  setwd(path)
  x<-list.files()
  x<-x[grep(".fas",x)]
  tab.tab<-NULL
  basepairs<-NULL
  for(i in 1:length(x)){
    
    fas<-read.dna(x[i],format="fasta")
    fas<-fas[match(tree.b[[1]],rownames(fas)),]
    length<-ncol(fas)
    fas<-fas2ms2(fas)
    basepairs<-c(basepairs,length-fas[[2]])
    fas<-ms.to.DNAbin(fas[[1]],bp.length=length-fas[[2]])
    
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
  if(moments==T){
    observed<-c(apply(tab.tab,2,mean),apply(tab.tab,2,var),apply(tab.tab,2,kurtosis),apply(tab.tab,2,skewness))
  } else {observed<-tab.tab}
  
  basepairs<-c(mean(basepairs),sd(basepairs))
  names(basepairs)<-c("basepairs mean and SD")
  return(list(t(data.frame(observed)),basepairs))
}




################### get spectral from a list of trees

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
## get a node in a phylogeny
get.node<-function(tree){
  output<-list()
  tree<-strsplit(tree,"")
  
  for(i in 1:length(tree[[1]])){
    if(tree[[1]][i]==")"){
      brack1<-i
      break}
  }
  for(i in brack1:1){
    if(tree[[1]][i]=="("){
      brack2<-i
      break}
  }
  
  for(i in brack1:1){
    if(tree[[1]][i]==","){
      comma<-i
      break}
  }
  
  # for(i in brack1:length(tree[[1]])){
  # if(tree[[1]][i]==","){
  #    comma2<-i
  #    break}
  #}
  
  junction<-tree[[1]][brack2:brack1]
  junction<-paste(junction,collapse="")
  t<-paste(tree[[1]],collapse="")
  output[[1]]<-gsub(junction,paste(tree[[1]][(comma+1):(brack1-1)],collapse=""),t,fixed=T)
  junction<-gsub("(","",junction,fixed=T)
  junction<-gsub(")","",junction,fixed=T)
  junction<-gsub(","," ",junction,fixed=T)
  #join.time<-paste(tree[[1]][(brack1+2):(comma2-1)],collapse="")
  output[[2]]<-junction
  #output[[3]]<-join.time
  return(output)
}
### get divergence time in the same order as nodes above
get.time<-function(tree,tree2){
  output<-list()
  
  node.matrix<-read.tree.nodes(paste(tree2[[1]],collapse=""))
  
  tree<-strsplit(tree,"")
  
  for(i in 1:length(tree[[1]])){
    
    if(tree[[1]][i]==")"){
      brack1<-i
      break}
  }
  for(i in brack1:1){
    if(tree[[1]][i]=="("){
      brack2<-i
      break}
  }
  
  for(i in brack1:1){
    if(tree[[1]][i]==","){
      comma<-i
      break}
  }
  
  for(i in comma:1){
    if(tree[[1]][i]==":"){
      comma2<-i
      break}
  }
  
  junction<-tree[[1]][brack2:brack1]
  junction<-paste(junction,collapse="")
  t<-paste(tree[[1]],collapse="")
  output[[1]]<-gsub(junction,paste(tree[[1]][(comma+1):(brack1-1)],collapse=""),t,fixed=T)
  junction<-gsub("(","",junction,fixed=T)
  junction<-gsub(")","",junction,fixed=T)
  nodes<-strsplit(junction,",")[[1]]
  junction<-gsub(","," ",junction,fixed=T)
  output[[2]]<-junction
  
  output[[3]]<-coaltime(which(node.matrix$names==getrid(nodes[1])),which(node.matrix$names==getrid(nodes[2])),
                        nodematrix = node.matrix$nodes,nspecies=length(node.matrix$names))
  
  return(output)
}
## get rid of branch lengths
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
## get three details
get.tree.info<-function(tree){
  
  tree.b<-getrid(tree)
  tree.b<-taxons.to.numbers(tree.b)
  
  ### get all joints of the tree
  input<-tree.b[[2]]
  join<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.node(input)
    join<-rbind(join,xx[[2]])
    input<-xx[[1]]
  }
  
  ### get all joint times of the tree
  input<-tree
  tree2<-tree
  join.t<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.time(input,tree2)
    join.t<-rbind(join.t,xx[[3]])
    input<-xx[[1]]
  }
  join<-cbind(join,join.t)
  tree.b[[3]]<-join
  return(tree.b)
}

## simualte network
sim.sp.tree<-function(tree,
                      Ne.prior,
                      biffurcating=T,
                      migration=F,
                      admixture=F,
                      mig,
                      hib.clade,
                      hib.priors,
                      major.sister,
                      minor.sister,
                      bp,
                      mi,
                      nsims,
                      ntrees,
                      segsites,
                      gen.time,
                      time.modif,
                      coaltrees,
                      spectral,
                      distance,
                      time.scalar=1)
{
  
  #### exclude branch lengths and change taxon to numbers
  tree.b<-getrid(tree)
  tree.b<-taxons.to.numbers(tree.b)
  
  ### get all joints of the tree
  input<-tree.b[[2]]
  join<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.node(input)
    join<-rbind(join,xx[[2]])
    input<-xx[[1]]
  }
  
  ### get all joint times of the tree
  input<-tree
  tree2<-tree
  join.t<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.time(input,tree2)
    join.t<-rbind(join.t,xx[[3]])
    input<-xx[[1]]
  }
  
  ####### generate ej (nodes) flag string. Nodes ages are rescaled.
  ej<-cbind(join.t,join)
  ej<-cbind(rep("-ej",nrow(ej)),ej)
  ej[,2]<-as.numeric(ej[,2])*time.scalar
  
  ###### generate en (ancestral Ne) flag. Sample random Nes. Ne values will change in the simulation cicle
  en<-cbind(rep("-en",nrow(ej)),ej[,2],ej[,3],runif(nrow(ej),Ne.prior[1],Ne.prior[2]))
  for(i in 1:nrow(en)){
    en[i,3]<-strsplit(ej[i,3]," ")[[1]][2]
  }
  
  #######
  ms.string<-list()
  simulated<-NULL
  ## master Ne
  master.Ne<-mean(Ne.prior)
  ## get nodes
  nodes<-as.numeric(ej[,2])
  
  #######################
  #######################
  # start simulations ###
  for(j in 1:nsims){
    sim.t<-NULL
    trees<-list()
    
    # put the original dates back on time strings
    ej[,2]<-nodes
    en[,2]<-nodes
    # rescale to coalescent (Ne proportion)
    ej[,2]<-as.numeric(ej[,2])/(4*master.Ne)/gen.time
    en[,2]<-as.numeric(ej[,2])#*1.01
    # sample time modifier
    time.mod<-runif(1,time.modif[1],time.modif[2])
    # modify time
    ej[,2]<-as.numeric(ej[,2])*time.mod
    en[,2]<-as.numeric(en[,2])*time.mod
    
    # sample an Ne mean
    Ne.mean<-runif(1,Ne.prior[1],Ne.prior[2])
    # sample Ne Standart deviation
    Ne.SD<-Ne.mean*runif(1,0.1,1)
    # sample contemporary Nes (truncated in 100 individuals)
    Nes<-rtnorm((nrow(ej)+1),Ne.mean,Ne.SD,lower=100)
    # Sample ancestral Nes (truncatedd 100 individuals)
    anc.Nes<-rtnorm(nrow(ej),Ne.mean,Ne.SD,lower=100)
    if(migration==T){
      Mig.prop<-runif(1,mig[1],mig[2])
      ntrees.mig<-floor(ntrees*Mig.prop)
      }else{
        Mig.prop<-0
      }
    if(admixture==T){
      Admix.prob.minor<-runif(1,hib.priors[4],hib.priors[5])
    }else{
      Admix.prob.minor<-0
    }
    ### simulate n-trees
    count<-0
    for(i in 1:ntrees){
      master.theta<-0
      while(master.theta<0.000001){
      # sample mutation rate
      mi.mean<-runif(1,mi[1],mi[2])
      mi.SD<-mi.mean*runif(1,0.1,1)
      #mi.SD<-mi.mean*mi.SD
      rate<-rnorm(1,mi.mean,mi.SD)
      while(rate<=0){
        rate<-rnorm(1,mi.mean,mi.SD)
      }
      # sample sequence length
      seq.length<-rnorm(1,bp[1],bp[2])
      # generate theta
      master.theta<-master.Ne*4*seq.length*(rate*gen.time)
      }
      # theta string
      ms.string[[1]]<-paste("-t",master.theta)
      ### pop structure
      ms.string[[2]]<-paste(c("-I",(nrow(ej)+1),(rep(1,nrow(ej)+1))),collapse=" ")
      ### current popsize strings
      n<-cbind(rep("-n",(nrow(ej)+1)),c(1:(nrow(ej)+1)),rep(0,(nrow(ej)+1)),Nes)
      n[,3]<-as.numeric(n[,4])/master.Ne
      ms.string[[3]]<-paste(apply(n[,1:3],1,paste,collapse=" "),collapse=" ")
      ### ancestral pop sizes string
      en[,4]<-anc.Nes
      en[,4]<-as.numeric(en[,4])/master.Ne
      ms.string[[4]]<-paste(apply(en,1,paste,collapse=" "),collapse=" ")
      ###  node string
      ms.string[[5]]<-paste(apply(ej,1,paste,collapse=" "),collapse=" ")
      
      
      #### network string
      
      if(biffurcating==F){
        if(admixture==T){
        split.time<-runif(1,hib.priors[1],hib.priors[2])*time.scalar
        ej.hib<-runif(1,split.time,(hib.priors[2]*time.scalar))
        split.time<-((split.time/gen.time)*time.mod)/(4*master.Ne)
        ej.hib<-((ej.hib/gen.time)*time.mod)/(4*master.Ne)
        es<-paste("-es",split.time,max(hib.clade),(1-Admix.prob.minor))
        ejh<-paste("-ej",ej.hib,(length(tree.b[[1]])+1),max(minor.sister))
        ms.string[[6]]<-paste(es,ejh)
      }
      if(migration==T){
        count<-count+1
        if(count>ntrees.mig){
          Mig<-0
        }else{
          Mig=1
        }
        #Mig<-rtnorm(1,Mig.mean,Mig.SD,0)
        mig.time<-((hib.priors[1]/gen.time)*time.mod)/(4*master.Ne)
        em<-paste("-em",mig.time,max(minor.sister),max(hib.clade),Mig)
        ms.string[[6]]<-em
      }
      }
      
      ######
      # combining all string peaces in one ms string
      ms.string.final<-paste(unlist(ms.string),collapse=" ")
      #print(ms.string.final)
      ########################################
      ########################################
      ########################################
      #### get the coalescent trees
      if(coaltrees==T){
        t<-ms(nreps = 1, nsam=(nrow(ej)+1),opts=paste("-T",ms.string.final))
        t<-strsplit(t,"//")[[3]]
        t<-read.tree(text=t)
        trees[[i]]<-t
      }
      
      #### build trees from simulated segregating sites
      if(segsites==T){
        fas<-ms.to.DNAbin(ms(nreps = 1, nsam=(nrow(ej)+1),opts=ms.string.final),bp.length = 0)
        
        while(length(fas)==0){
          fas<-ms.to.DNAbin(ms(nreps = 1, nsam=(nrow(ej)+1),opts=ms.string.final),bp.length = 0)  
        }
        
        d<-dist.dna(fas, model="N")/seq.length
        
        if(distance==T){
          sim.t<-rbind(sim.t,as.vector(d))
        }
       
        if(spectral==T){
          upgma.tree<-upgma(d)
          upgma.tree$edge.length<-upgma.tree$edge.length*1000
          trees[[i]]<-upgma.tree
        }
        
      }
      rm(fas)
    }
    
    if(spectral==T){
      eigen<-get.panda(trees,method="standard")
      sim.t<-cbind(eigen$principal_eigenvalue,eigen$asymmetry,eigen$eigenvalues)
      if(nrow(sim.t)>1){
        sim<-c(apply(sim.t,2,mean),apply(sim.t,2,var),apply(sim.t,2,kurtosis),apply(sim.t,2,skewness))
      } else {sim<-t(sim.t)} 
    }
    if(distance==T){
      if(segsites==F){
        stop("you need to have segsites=TRUE to get the distance")} else {
          if(nrow(sim.t)>1){
            nam<-t(combn(attr(d,"Labels"),2))
            nam<-apply(nam,1,paste,collapse="_")
            colnames(sim.t)<-nam
            sim<-c(apply(sim.t,2,mean))#,apply(sim.t,2,var),apply(sim.t,2,kurtosis),apply(sim.t,2,skewness))
            mean<-paste("mean",names(sim[1:length(nam)]),sep="_")
            #var<-paste("var",names(sim[(length(nam)+1):(2*length(nam))]),sep="_")
            #kur<-paste("kur",names(sim[((2*length(nam))+1):(3*length(nam))]),sep="_")
            #skew<-paste("skew",names(sim[((3*length(nam))+1):(4*length(nam))]),sep="_")
            names(sim)<-c(mean)#,var,kur,skew)
          } else {
            nam<-t(combn(attr(d,"Labels"),2))
            nam<-apply(nam,1,paste,collapse="_")
            colnames(sim.t)<-nam
            sim<-t(sim.t)} 
        }
    }
    #print(i)
    sim<-cbind(time.mod,Ne.mean,Ne.SD,mi.mean,mi.SD,Mig.prop,Admix.prob.minor,t(sim))
    simulated<-rbind(simulated,sim)
    rm(time.mod,Ne.mean,Ne.SD,mi.mean,mi.SD,Mig.prop,Admix.prob.minor,sim)
    print(j)
  }
  return(simulated)
}

