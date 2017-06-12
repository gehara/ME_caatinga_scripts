
auto.BSLP<-function(path.to.fastas){
setwd(path.to.fastas)
files<-list.files()
files<-files[c(grep(".fas",files))]
files<-c(files,files[c(grep(".fasta",files))])
for(i in 1:length(files)){
  x<-strsplit(files[i],".fas")[[1]]
  system(paste("java -jar beastgen.jar BSP.template ",files[i]," > ",x,".xml",sep=""))
  command<-paste("java -jar beast.jar -overwrite -beagle_off ",x,".xml",sep="")
  system(command)
}
}
auto.BSLP(getwd())