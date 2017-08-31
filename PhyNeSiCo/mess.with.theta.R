
## this function put thetas in a initial input, a nwk string with branch lengths.
## all thetas are places after branch lengths and they all equal 0.01 at that stage.
## the nwk string needs to have branch length information!

put.thetas.on.nwk<-function(tree){
  x<-strsplit(tree,":")[[1]]
  for(i in 2:length(x)){
    y<-strsplit(x[i],"")[[1]]
    for(j in 1:length(y)){
      if(y[j] %in% c(",",")",";")){
        x[i]<-paste(paste(y[1:(j-1)],collapse=""),"#0.01",paste(y[j:length(y)],collapse=""),sep="")
        next
      } 
    }
    }
  x<-paste(x,collapse=":")
  return(x)
}


### this is part of the main function bellow and
# is used to get the thetas from the tree string 
# ex of string input: "(((S1:0.005#0.01,S2:0.005#0.01):0.00025#.01,S3:0.00525#0.01):0.00025#0.01,S4:0.0055#0.01)#.01;"

get.thetas<-function(tree){
  tree<-strsplit(tree,"#")[[1]]
  thetas<-NULL
  for(i in 2:length(tree)){
    
    x<-strsplit(tree[i],"")[[1]]
    
    for(j in 1:length(x)){
      
      if(x[j] %in% c(",",")",";")){
        t<-paste(x[1:(j-1)],collapse="")
      }
    }
    thetas<-c(thetas,t)
  }
  return(thetas)
}

# this function replace thetas with new sampled values from a gamma distribution
# It uses the function above to get the thetas from the tree string.
# You need a tree string as input and two gamma parameters (gamma.shape and gamma.rate)

# It returns a list with the new tree string as the first element, and the new thetas as a second element.

mess.with.thetas<-function(tree,gamma.shape,gamma.rate){
  
  thetas<-get.thetas(tree)
  new.thetas<-round(rgamma(length(thetas),gamma.shape,gamma.rate),4)
  tree<-strsplit(tree,"#")[[1]]
  for(i in 1:length(thetas)){
    tree[i+1]<-gsub(thetas[i],new.thetas[i],tree[i+1])
    }
  tree<-paste(tree,collapse="#")
  output<-list()
  output$tree<-tree
  output$thetas<-new.thetas
  return(output)
}


### Example!
#string without theta. only branch length info
tree<-"(Cemophora_coccinea:16.3999643,((Lampropeltis_extenuata:5.466656054,((Lampropeltis_rhombomaculata:1.822230181,Lampropeltis_calligaster:1.822230181):1.822200026,((Lampropeltis_getula:1.214816654,Lampropeltis_nigra:1.214816654):1.214823732,(Lampropeltis_holbrooki:1.214519523,(Lampropeltis_californiae:0.6071318993,Lampropeltis_splendida:0.6071318993):0.6073876235):1.215120863):1.214789821):1.822225848):5.466689376,(Lampropeltis_elapsoides:5.46665917,((Lampropeltis_annulata:1.822219939,Lampropeltis_gentilis:1.822219939):1.822217324,((Lampropeltis_polyzona:1.214541049,(Lampropeltis_abnorma:0.6071502151,Lampropeltis_micropholis:0.6071502151):0.6073908336):1.214935994,((Lampropeltis_webbi:0.8094465841,(Lampropeltis_mexicana:0.4045359888,Lampropeltis_ruthveni:0.4045359888):0.4049105953):0.8101756765,(Lampropeltis_alterna:0.8093582685,((Lampropeltis_zonata:0.2694659048,Lampropeltis_multifasciata:0.2694659048):0.2699155957,(Lampropeltis_knoblochi:0.2694549107,Lampropeltis_pyromelana:0.2694549107):0.2699265898):0.269976768):0.8102639922):0.8098547822):1.21496022):1.822221907):5.46668626):5.466618865);"

tree<-put.thetas.on.nwk(tree)

mess.with.thetas(tree,gamma.shape = 1,gamma.rate = 10)

