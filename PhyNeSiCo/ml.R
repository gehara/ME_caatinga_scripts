library(gbm)
library(caret)
library(pROC)
library(randomForest)

setwd("~/Dropbox/AMNH/Lampropeltis/Lampropeltis_TOL/Code")
source("network.sim.R")

bif.beast<-read.table("bif_beast_distance.txt")
net.beast<-read.table("net_beast_distance.txt")
#bif.astral<-read.table("bif_astral_distance.txt")
#net.astral<-read.table("net_astral_distance.txt")
head(bif.beast)
#colnames(x)->colnames(bif.astral)
colnames(x)->colnames(bif.beast)
#colnames(x)->colnames(net.astral)
colnames(x)->colnames(net.beast)

setwd("/home/marcelo/Dropbox/AMNH/Lampropeltis/gene_fasta_11taxa/")
astral.tree<-"(((((zonata:0.2694659048,multifasciata:0.2694659048):0.2699155957,(knoblochi:0.2694549107,pyromelana:0.2694549107):0.2699265898):0.269976768,alterna:0.8093582685):0.8102639922,((mexicana:0.4045359888,ruthveni:0.4045359888):0.4049105953,webbi:0.8094465841):0.8101756765):0.8098547822,((abnorma:0.6071502151,micropholis:0.6071502151):0.6073908336,polyzona:1.214541049):1.214935994);"
observed.dist<-get.observed(tree.b=get.tree.info(astral.tree),path=getwd(),spectral=F, distance=T, momentum=T)
observed.dist<-t(data.frame(observed.dist))

observed.dist<-observed.dist[,1:55]

index<-c(rep("bif",nrow(bif.beast)),rep("net",nrow(net.beast)))
      
#models.astral<-rbind(bif.astral[,6:ncol(bif.astral)],net.astral[,6:ncol(net.astral)])

models.beast<-rbind(bif.beast[,6:60],net.beast[,6:60])
models.beast<-data.frame(models.beast)

colnames(models.beast)->colnames(observed.dist)
models.beast<-data.frame(models.beast)
head(models.beast)

descrCor <- cor(models.beast)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .90)
models.uncor <- models.beast[,-highlyCorDescr]
head(models.uncor)

head(observed.dist)
observed.uncor<-observed.dist[-highlyCorDescr]
t(data.frame(observed.uncor))->observed.uncor
colnames(models.uncor)->colnames(observed.uncor)

head(models.uncor)
head(observed.uncor)

#descrCor2 <- cor(models)
#summary(descrCor2[upper.tri(descrCor2)])


delta<-NULL
for(i in 1:50){

  #models <- models.astral[,-highlyCorDescr]
  models.beast<-cbind(models.beast,index)

  outcomeName <- 'index'
  predictorsNames <- names(models.beast)[names(models.beast) != outcomeName]

  splitIndex <- createDataPartition(models.beast[,outcomeName], p = .75, list = FALSE, times = 1)
  train <- models.beast[ splitIndex,]
  test  <- models.beast[-splitIndex,]

  objControl <- trainControl(method='boot', number=10, returnResamp='final', 
                           summaryFunction = twoClassSummary, classProbs = TRUE)
  
  nnetModel <-train(train[,predictorsNames], train[,outcomeName], 
                  method="nnet",maxit=1000, 
                  trControl=objControl,  
                  metric = "ROC",
                  preProc = c("center", "scale"))

  predictions <- predict(object=nnetModel, test[,predictorsNames], type='raw')
  accu<-postResample(pred=predictions, obs=as.factor(test[,outcomeName]))[1]

    #zz<-test[sample(c(1:nrow(test)),nrow(test)),1:ncol(test)-1]
    #test<-cbind(zz,test[,ncol(test)])
    #test<-data.frame(test)
    #colnames(test)<-c(colnames(test)[1:28],"index")
    #predictions <- predict(object=nnetModel, test[,predictorsNames], type='raw')
    #accu.random<-postResample(pred=predictions, obs=as.factor(test[,outcomeName]))[1]

    #d<-accu-accu.random
    #delta<-c(delta,d)
}
head(observed.dist)
head(models.beast)

predict(object=nnetModel, t(data.frame(observed.dist[,1:55])), type='prob')


getwd()

write.table(delta,"delta.txt", quote=F)


nnetpcaModel <- train(train[,predictorsNames], train[,outcomeName], 
                      method="nnet", maxit=1000
                      trControl=objControl,  
                      metric = "ROC",
                      preProc = c("center", "scale","pca"))

rfModel <- train(train[,predictorsNames], train[,outcomeName], 
                 method="rf",maxit=1000,
                 trControl=objControl,  
                 metric = "ROC",
                 preProc = c("center", "scale"))

observed<-t(data.frame(observed.dist[,1:55]))

prediction<-predict(object=nnetModel, observed, type='prob')
prediction<-rbind(prediction,predict(object=nnetpcaModel, observed, type='prob'))
prediction<-rbind(prediction,predict(object=rfModel, observed, type='prob'))


predictions <- predict(object=nnetModel, test[,predictorsNames], type='raw')
print(postResample(pred=predictions, obs=as.factor(test[,outcomeName])))

predictions <- predict(object=nnetModel, test[,predictorsNames], type='prob')
auc <- roc(ifelse(test[,outcomeName]=="net",1,0), predictions[[2]])
AUC <- auc$auc

predictions <- predict(object=nnetpcaModel, test[,predictorsNames], type='prob')
auc <- roc(ifelse(test[,outcomeName]=="net",1,0), predictions[[2]])
AUC <- rbind(AUC,auc$auc)

predictions <- predict(object=rfModel, test[,predictorsNames], type='prob')
auc <- roc(ifelse(test[,outcomeName]=="net",1,0), predictions[[2]])
AUC <- rbind(AUC,auc$auc)

prediction<-cbind(prediction,AUC)

rownames(prediction)<-c("NNET","NNETpca","RF")

write.table(prediction,"prediction_astral.txt", quote=F, row.names = T, col.names = T, sep="\t")



summary(nnetModel)
Imp <- varImp(nnetModel, scale = FALSE)
plot(gbmImp, top = 20)
plot.nnet(nnetModel)

predictions <- predict(object=nnetpcaModel, test[,predictorsNames], type='raw')
print(postResample(pred=predictions, obs=as.factor(test[,outcomeName])))

predictions <- predict(object=nnetModel, test[,predictorsNames], type='prob')
auc <- roc(ifelse(test[,outcomeName]=="net",1,0), predictions[[2]])
auc$auc
predict(object=nnetpcaModel, observed, type='prob')



AUC<-NULL
for(i in 1:100){
  models$index<-sample(models$index,nrow(models))
  
  train <- models[ splitIndex,]
  test  <- models[-splitIndex,]
  
  nnetModel <- train(train[,predictorsNames], train[,outcomeName], 
                     method="nnet",maxit=1000,
                     trControl=objControl,  
                     metric = "ROC",
                     preProc = c("center", "scale"))
  
  predictions <- predict(object=nnetModel, test[,predictorsNames], type='prob')
  auc <- roc(ifelse(test[,outcomeName]=="bif",1,0), predictions[[2]])
  AUC<-c(AUC,auc$auc)
  print(i)
}

plot(density(main.AUC-AUC))

nnetpcaModel <- train(train[,predictorsNames], train[,outcomeName], 
                  method="nnet", 
                  trControl=objControl,  
                  metric = "ROC",
                  preProc = c("center", "scale","pca"))

rfModel <- train(train[,predictorsNames], train[,outcomeName], 
                   method="rf",maxit=1000,
                   trControl=objControl,  
                   metric = "ROC",
                   preProc = c("center", "scale"))

rfpcaModel <- train(train[,predictorsNames], train[,outcomeName], 
                      method="rf", 
                      trControl=objControl,  
                      metric = "ROC",
                      preProc = c("center", "scale","pca"))



summary(objModel)
gbmImp <- varImp(objModel, scale = FALSE)
plot(gbmImp, top = 20)
print(objModel)
?nnet

predictions <- predict(object=nnetModel, test[,predictorsNames], type='raw')
print(postResample(pred=predictions, obs=as.factor(test[,outcomeName])))
predictions <- predict(object=nnetModel, test[,predictorsNames], type='prob')
head(predictions)
auc <- roc(ifelse(test[,outcomeName]=="bif",1,0), predictions[[2]])
print(auc$auc)





colnames(observed.dist)<-colnames(models[,1:220])
observed.dist[,-highlyCorDescr]-> observed

ncol(models)
ncol(observed)



predict(object=nnetModel, observed, type='prob')
predict(object=nnetpcaModel, observed, type='prob')
predict(object=rfModel, observed, type='prob')

#predict(object=rfpcapcaModel, observed, type='prob')

predictions <- predict(object=objModel, observed, type='prob')
predictions

getModelInfo()$glmnet$type


plot(density(auc$auc-accu))
