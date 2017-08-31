#example showing use of separate colors for input layer
#color based on relative importance using 'gar.fun'

##
#create input data
seed.val<-3
set.seed(seed.val)

num.vars<-8
num.obs<-1000

#input variables
library(clusterGeneration)
cov.mat<-genPositiveDefMat(num.vars,covMethod=c("unifcorrmat"))$Sigma
rand.vars<-mvrnorm(num.obs,rep(0,num.vars),Sigma=cov.mat)

#output variables
parms<-runif(num.vars,-10,10)
y1<-rand.vars %*% matrix(parms) + rnorm(num.obs,sd=20)

#final datasets
rand.vars<-data.frame(rand.vars)
resp<-data.frame(y1)
names(resp)<-'Y1'
dat.in<-data.frame(resp,rand.vars)

##
#create model
library(nnet)
mod1<-nnet(rand.vars,resp,data=dat.in,size=10,linout=T)

##
#relative importance function
library(devtools)
library(reshape)
install.packages("reshape")
source_url('https://gist.github.com/fawda123/6206737/raw/2e1bc9cbc48d1a56d2a79dd1d33f414213f5f1b1/gar_fun.r')

#relative importance of input variables for Y1
rel.imp<-gar.fun(predictorsNames,nnetModel$finalModel,bar.plot=F)$rel.imp

nnetModel$results

#color vector based on relative importance of input values
cols<-colorRampPalette(c('green','red'))(69)[rank(rel.imp)]

##
#plotting function
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')

#plot model with new color vector
#separate colors for input vectors using a list for 'circle.col'
plot(nnetModel$finalModel,circle.col=list(cols,'lightblue'))
