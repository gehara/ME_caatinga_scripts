


setwd("~/Dropbox/AMNH/Lampropeltis/TICR/my-genes-mb/trees.nwk")
library(ape)
x<-list.files()
tre<-read.tree(x[307])

library(pegas)
library(phytools)
library(geiger)

group<-tre[[1]]$tip.label
group<-group[c(14,16,11,19,20,17,18,15,12,13,22)]

group1<-group[c(1:8)]
group2<-group[c(6:11)]

A<-NULL
B<-NULL
for(i in 1:length(tre)){
  a<-gsi(tre[[i]],group1)
  b<-gsi(tre[[i]],group2)
  A<-c(A,a)
  B<-c(B,b)
  print(i)
}
table<-cbind(A,B)
y<-NULL
for(i in 1:nrow(table)){
  x<-max(table[i,])
  y<-c(y,which(table[i,]==x))
}
length(which(y==1))/nrow(table)
length(which(y==2))/nrow(table)
#####################
w<-NULL
z<-NULL  
for(j in 1:100){
A<-NULL
B<-NULL
for(i in 1:length(tre)){
  group.random<-match(group,tre[[i]]$tip.label)
  tre[[i]]$tip.label[group.random]<-sample(tre[[i]]$tip.label[group.random])
  
  a<-gsi(tre[[i]],group1)
  b<-gsi(tre[[i]],group2)
  A<-c(A,a)
  B<-c(B,b)
  #print(i)
}
table<-cbind(A,B)
y<-NULL
for(i in 1:nrow(table)){
  x<-max(table[i,])
  y<-c(y,which(table[i,]==x))
}
z<-c(z,length(which(y==1))/nrow(table))
w<-c(w,length(which(y==2))/nrow(table))
print(j)
}

plot(density(w))
plot(density(z))
length(which(y==1))/nrow(table)
length(which(y==2))/nrow(table)

