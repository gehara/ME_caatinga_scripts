source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")

treeA<-"(Cemophora_coccinea:11.678229466235559,(((Lampropeltis_elapsoides:6.240913466415052,(Lampropeltis_annulata:2.1388501466619054,Lampropeltis_gentilis:2.1388501466619054):4.102063319753147):2.3046166642738006,(((Lampropeltis_abnorma:2.413413791353035,Lampropeltis_micropholis:2.413413791353035):2.046516169196659,Lampropeltis_polyzona:4.459929960549694):2.5425428142366115,((Lampropeltis_alterna:4.554400812401852,((Lampropeltis_zonata:1.5346400910659765,Lampropeltis_multifasciata:1.5346400910659765):1.4232385476943845,(Lampropeltis_knoblochi:1.6340343424803274,Lampropeltis_pyromelana:1.6340343424803274):1.3238442962800336):1.5965221736414925):1.2448383594795835,(Lampropeltis_webbi:2.892264576497898,(Lampropeltis_mexicana:1.539058444509318,Lampropeltis_ruthveni:1.539058444509318):1.3532061319885798):2.9069745953835375):1.20323360290487):1.5430573559025476):1.9194000035200798,(((Lampropeltis_nigra:2.623141453795378,Lampropeltis_getula:2.623141453795378):1.7038376207295447,(Lampropeltis_holbrooki:2.969604917867215,(Lampropeltis_californiae:1.6295385079701354,Lampropeltis_splendida:1.6295385079701354):1.3400664098970814):1.3573741566577078):2.3942359444784973,(Lampropeltis_extenuata:5.480392736726901,(Lampropeltis_rhombomaculata:2.3301760341467883,Lampropeltis_calligaster:2.3301760341467883):3.1502167025801127):1.2408222822765191):3.7437151152055126):1.2132993320266259);"
treeA<-read.tree(text=treeA)
treeB<-"(Cemophora_coccinea:11.646531622834841,(((Lampropeltis_elapsoides:6.4724541318431905,(Lampropeltis_annulata:2.223557678670028,Lampropeltis_gentilis:2.223557678670028):4.248896453173162):1.7351824579374346,(((Lampropeltis_abnorma:2.299446356630222,Lampropeltis_micropholis:2.299446356630222):3.015957032048095,(Lampropeltis_polyzona:4.494348796582882,(Lampropeltis_webbi:2.437111922702547,(Lampropeltis_mexicana:1.2095584711996956,Lampropeltis_ruthveni:1.2095584711996956):1.2275534515028532):2.057236873880335):0.8210545920954351):1.5441590364726077,(Lampropeltis_alterna:4.98585614680508,((Lampropeltis_zonata:1.6551188768199818,Lampropeltis_multifasciata:1.6551188768199818):1.498265321401723,(Lampropeltis_knoblochi:1.7104033936538254,Lampropeltis_pyromelana:1.7104033936538254):1.4429808045678794):1.8324719485833754):1.8737062783458445):1.3480741646297005):2.0683784200558453,(((Lampropeltis_nigra:2.7940973017001376,Lampropeltis_getula:2.7940973017001376):1.6979135085595853,(Lampropeltis_holbrooki:3.1538217772605943,(Lampropeltis_californiae:1.7761760182728263,Lampropeltis_splendida:1.7761760182728263):1.3776457589877698):1.3381890329991286):2.4331913359368373,(Lampropeltis_extenuata:5.663054476446628,(Lampropeltis_rhombomaculata:2.453455947146429,Lampropeltis_calligaster:2.453455947146429):3.209598529300199):1.2621476697499325):3.35081286363991):1.3705166129983706);"  
treeB<-read.tree(text=treeB)
plot(treeA)
plot(treeB)

#treeA$edge.length<-treeA$edge.length/1000
#treeB$edge.length<-treeB$edge.length/1000

for(i in 1:length(treeA$tip.label)){
  treeA$tip.label[i]<-strsplit(treeA$tip.label,"_")[[i]][2]
}

for(i in 1:length(treeB$tip.label)){
  treeB$tip.label[i]<-strsplit(treeB$tip.label,"_")[[i]][2]
}

out<-treeA$tip.label[-match(group,treeA$tip.label)]

sub.treeA<-drop.tip(treeA,out)
sub.treeB<-drop.tip(treeB,out)
plot(sub.treeA)
plot(sub.treeB)

sub.treeA$edge.length<-sub.treeA$edge.length/1000000
tA<-write.tree(sub.treeA)
plot(sub.treeA)

tA<-getrid(tA)
tA<-taxons.to.numbers(tA)


tA.matrix<-read.tree.nodes(tA)

species.name(tA)

tA.matrix$nodes[order(tA.matrix$nodes[,4]),]

J.matrix<-cbind(c(1:11),tA.matrix$nodes[1:11,])

match()

N.nodes<-length(unique(tA.matrix$nodes[,1]))

nodes<-unique(tA.matrix$nodes[,1])


all.join<-NULL
for(i in 1:length(nodes)){
  join<-c(which(tA.matrix$nodes[,1]==nodes[i]),nodes[i])
  all.join<-rbind(all.join,join)
}

J.matrix<-tA.matrix$nodes[12:nrow(tA.matrix$nodes),]

J.matrix[order(J.matrix[,4]),]
