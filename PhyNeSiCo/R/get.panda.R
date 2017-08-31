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