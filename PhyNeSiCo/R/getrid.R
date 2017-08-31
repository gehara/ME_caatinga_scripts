## get rid of branch lengths
getrid<-function(x){
  z<-c(0:9,".",":")
  for(i in z){
    x<-gsub(i,"",x,fixed=T)
  }
  return(x)
}