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