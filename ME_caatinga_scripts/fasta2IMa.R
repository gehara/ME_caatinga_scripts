
fasta2IMa2<-function(path,mi){
library(ape)
setwd(path)
files<-list.files()
fasta.files<-files[grep(".fas",files)]

for(u in 1:length(fasta.files)){
  fas<-read.dna(file=fasta.files[u], format="fasta")
  seqnames<-rownames(fas)
  fas<-as.character(fas)
  seqlength<-ncol(fas)
  nsamp<-nrow(fas)
  sequence<-vector()
  for(i in 1:nrow(fas)){
    sequence[i]<-paste(paste("sample"),paste(fas[i,],collapse=""),sep="    ")
  }
  write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".IMa2",sep=""),"#IMa input")
  write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".IMa2",sep=""),"1", append=T)
  write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".IMa2",sep=""),"pop1", append=T)
  write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".IMa2",sep=""),"1", append=T)
  write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".IMa2",sep=""),"1", append=T)
  write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".IMa2",sep=""),paste("locus1",nsamp,seqlength,"H",1,(seqlength*mi)),append=T)
  write(file=paste(strsplit(fasta.files[u],".fas")[[1]],".IMa2",sep=""),paste(sequence,sep="\n"),append=T)
}
}
fasta2IMa2(mi=1e-8,path=getwd())





