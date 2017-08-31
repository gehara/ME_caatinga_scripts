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