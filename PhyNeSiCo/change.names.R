setwd("/home/marcelo/Dropbox/AMNH/Lampropeltis/TICR/my-genes-mb/trees.nwk/")
library(ape)
x<-read.tree("list")
for(j in 1:length(x)){
y<-x[j][[1]]$tip.label
y<-gsub("Cemorphora","Cemophora",y)
y<-gsub("calligaster_east","rhombomaculata",y)
y<-gsub("calligaster_west","calligaster",y)
y<-gsub("triangulum_Honduras","abnorma",y)
y<-gsub("triangulum_Mex","polyzona",y)
y<-gsub("triangulum_TAM","annulata",y)
y<-gsub("triangulum_SA","micropholis",y)
y<-gsub("triangulum_west","gentilis",y)
y<-gsub("zonata_south","multifasciata",y)
y<-gsub("zonata_north","zonata",y)
y<-gsub("niger","nigra",y)
y<-gsub("extenuatum","extenuata",y)
for(i in 1:length(y)){
  y[i]<-paste(strsplit(y[i],"_",fixed=T)[[1]][3:4],collapse = "_")
}
x[j][[1]]$tip.label<-y
x[j][[1]]$edge.length<-NULL
}
for(i in 1:length(x)){
x[[i]]<-drop.tip(x[[i]],tip=y[c(1:3,7:13,20:21)])
}

write.tree(x,file="subclade1_all_trees")
library(help=ape)

x<-readLines(con="subclade1_all_trees")

for(i in 1:length(x)){
x[i]<-paste("Tree ","gene",i," =  ",x[i],sep="")
}
write(x,file="subclade1_all_trees")

