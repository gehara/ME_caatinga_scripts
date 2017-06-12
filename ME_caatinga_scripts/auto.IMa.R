auto.IMa2<-function(path=getwd()){
  library(ggplot2)
  setwd(path)
  files<-list.files()
  IMa.files<-files[grep(".IMa2",files)]

  for(i in 1:length(IMa.files)){
    x<-strsplit(IMa.files[i],".IMa2")[[1]]
    m.rate<-readLines(IMa.files[i])[6]
    m.rate<-as.numeric(strsplit(m.rate," ")[[1]][6])
    system(paste("./IMa2 -i",IMa.files[i]," -o",x,".out -p23 -q50 -b200000 -l10000",sep=""))
  
    out<-readLines(paste(x,".out",sep=""))
    line1<-grep("q0 curve",out)
    line2<-grep("0_Ka curve",out)-1
    qcurve<-out[line1:line2]
    write(qcurve,file=paste(x,".qcurve.txt",sep=""))
  
    line1<-grep("Parameter	q0	P",out,fixed=T)+2
    data<-read.table(paste(x,".out",sep=""), skip=line1,nrows=1000)
    colnames(data)<-c("HiPt","Ne","p")
    data[,2]<-round((((data[,2]/m.rate)/4)/0.25/1000),0)
    row<-grep(max(data[,3]),data[,3])[[1]]
      f<-qplot(Ne, p, data=data, xlab="Thousend Individuals", ylab="p", main="(Ne x gentime) posterior distribution")+
        theme(text = element_text(size=10))+
        scale_x_continuous(breaks=seq(0,max(data[,2]),data[row,2]))
    pdf(file=paste(x,".Ne.pdf",sep=""),paper="a4r")
    plot(f)
    dev.off()
  
    line1<-grep("Parameter	g0_tmrca",out,fixed=T)+2
    data<-read.table(paste(x,".out",sep=""), skip=line1,nrows=1000)
    colnames(data)<-c("HiPt","t","p")
    row<-grep(max(data[,3]),data[,3])[[1]]
      f<-qplot(t, p, data=data, xlab="time", ylab="p", main="tau posterior distribution")+
          theme(text = element_text(size=10))+
          scale_x_continuous(breaks=seq(0,max(data[,2]),data[row,2]))
    pdf(file=paste(x,".tau.pdf",sep=""),paper="a4r")
    plot(f)
    dev.off()
  
    sum<-grep("Summaries",out,fixed=T)
    result<-list()
    result[[1]]<-read.table(paste(x,".out",sep=""), skip=sum[1],nrows=8, header=T)
    result[[2]]<-read.table(paste(x,".out",sep=""), skip=sum[2],nrows=9, header=T)
  
    Ne<-((result[[1]][,2]/m.rate)/4)/0.25
    result[[1]]<-cbind(result[[1]],Ne)
    colnames(result[[1]])<-c("Value","theta","Ne x gentime")
    tmrca<-result[[2]][,2]/m.rate
    result[[2]]<-cbind(result[[2]],tmrca)
    write.table(result[[1]],file=paste(x,".NeTable.txt",sep=""),quote=F,sep="\t",col.names = T,row.names = F)
    write.table(result[[2]],file=paste(x,".tmrcaTable.txt",sep=""),quote=F,sep="\t",col.names = T,row.names = F)
}
}
auto.IMa2()
