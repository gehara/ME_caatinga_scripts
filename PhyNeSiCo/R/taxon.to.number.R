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