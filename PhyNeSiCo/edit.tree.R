### read in tree major topology
### beast tree
treeA<-"(Cemophora_coccinea:11.678229466235559,(((Lampropeltis_elapsoides:6.240913466415052,(Lampropeltis_annulata:2.1388501466619054,Lampropeltis_gentilis:2.1388501466619054):4.102063319753147):2.3046166642738006,(((Lampropeltis_abnorma:2.413413791353035,Lampropeltis_micropholis:2.413413791353035):2.046516169196659,Lampropeltis_polyzona:4.459929960549694):2.5425428142366115,((Lampropeltis_alterna:4.554400812401852,((Lampropeltis_zonata:1.5346400910659765,Lampropeltis_multifasciata:1.5346400910659765):1.4232385476943845,(Lampropeltis_knoblochi:1.6340343424803274,Lampropeltis_pyromelana:1.6340343424803274):1.3238442962800336):1.5965221736414925):1.2448383594795835,(Lampropeltis_webbi:2.892264576497898,(Lampropeltis_mexicana:1.539058444509318,Lampropeltis_ruthveni:1.539058444509318):1.3532061319885798):2.9069745953835375):1.20323360290487):1.5430573559025476):1.9194000035200798,(((Lampropeltis_nigra:2.623141453795378,Lampropeltis_getula:2.623141453795378):1.7038376207295447,(Lampropeltis_holbrooki:2.969604917867215,(Lampropeltis_californiae:1.6295385079701354,Lampropeltis_splendida:1.6295385079701354):1.3400664098970814):1.3573741566577078):2.3942359444784973,(Lampropeltis_extenuata:5.480392736726901,(Lampropeltis_rhombomaculata:2.3301760341467883,Lampropeltis_calligaster:2.3301760341467883):3.1502167025801127):1.2408222822765191):3.7437151152055126):1.2132993320266259);"
### astral tree
treeA<-"(Cemophora_coccinea:16.3999643,((Lampropeltis_extenuata:5.466656054,((Lampropeltis_rhombomaculata:1.822230181,Lampropeltis_calligaster:1.822230181):1.822200026,((Lampropeltis_getula:1.214816654,Lampropeltis_nigra:1.214816654):1.214823732,(Lampropeltis_holbrooki:1.214519523,(Lampropeltis_californiae:0.6071318993,Lampropeltis_splendida:0.6071318993):0.6073876235):1.215120863):1.214789821):1.822225848):5.466689376,(Lampropeltis_elapsoides:5.46665917,((Lampropeltis_annulata:1.822219939,Lampropeltis_gentilis:1.822219939):1.822217324,((Lampropeltis_polyzona:1.214541049,(Lampropeltis_abnorma:0.6071502151,Lampropeltis_micropholis:0.6071502151):0.6073908336):1.214935994,((Lampropeltis_webbi:0.8094465841,(Lampropeltis_mexicana:0.4045359888,Lampropeltis_ruthveni:0.4045359888):0.4049105953):0.8101756765,(Lampropeltis_alterna:0.8093582685,((Lampropeltis_zonata:0.2694659048,Lampropeltis_multifasciata:0.2694659048):0.2699155957,(Lampropeltis_knoblochi:0.2694549107,Lampropeltis_pyromelana:0.2694549107):0.2699265898):0.269976768):0.8102639922):0.8098547822):1.21496022):1.822221907):5.46668626):5.466618865);"

treeA<-read.tree(text=treeA)
plot(treeA)
#### change names, exclude genus
for(i in 1:length(treeA$tip.label)){
  treeA$tip.label[i]<-strsplit(treeA$tip.label,"_")[[i]][2]
}

### subset the main hybrid clade
group<-c("alterna","zonata","multifasciata","pyromelana","knoblochi","webbi",
         "ruthveni","mexicana","polyzona","abnorma","micropholis")
out<-treeA$tip.label[-match(group,treeA$tip.label)]
treeA<-drop.tip(treeA,out)
plot(treeA)
write.tree(treeA)
