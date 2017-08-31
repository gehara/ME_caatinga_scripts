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