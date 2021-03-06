## importation of the main libraries
library(mixer)
library(ergm)
library(network)
library(igraph)
library("RColorBrewer") ### pour colorer nos graphes
library(flexclust)
library(lda)
library(sand)

setwd(dir="/Users/yoanrussac/desktop/ENSAE-Cours/2A/2nd Semestre/Séminaire de stat/StochasticBlockModel/Data")
dir()
### importation des données et traitement des bases de données
donnee_MMB=read.csv2(file='MMB_data.csv',sep=",",dec=".")
donnee_MMB=as.data.frame(donnee_MMB)
donnee_twitter=read.csv2(file='donnes_twitter.csv',sep=";",dec=".")
donnee_twitter$node<-as.factor(donnee_twitter$node)

donnee_twitter_1000=read.csv2(file='donnees_twitter_1000.csv',sep=";",dec=".")
donnee_twitter_1000$node<-as.factor(donnee_twitter_1000$node)




adja_donnee_twitter=read.csv2(file='adja_donnees_twitter.csv',sep=";",dec=".")
adja_donnee_twitter[,1]<-as.factor(adja_donnee_twitter[,1])
rownames(adja_donnee_twitter)=adja_donnee_twitter[,1]
adja_donnee_twitter=adja_donnee_twitter[,-1]
colnames(adja_donnee_twitter)=rownames(adja_donnee_twitter)

adja_donnee_twitter_1000=read.csv(file='adja_1000.csv',sep=";",dec=".")
rownames(adja_donnee_twitter_1000)=adja_donnee_twitter_1000[,1]
adja_donnee_twitter_1000=adja_donnee_twitter_1000[,-1]
colnames(adja_donnee_twitter_1000)=rownames(adja_donnee_twitter_1000)


indice_ordonne=c()
for (i in 1:length(rownames(adja_donnee_twitter))){
  print(which(donnee_twitter$node==rownames(adja_donnee_twitter)[i]))
  indice_ordonne=c(indice_ordonne,which(donnee_twitter$node==rownames(adja_donnee_twitter)[i]))
}

indice_ordonne_1000=c()
for (i in 1:length(rownames(adja_donnee_twitter_1000))){
  indice_ordonne_1000=c(indice_ordonne_1000,which(donnee_twitter_1000$user==rownames(adja_donnee_twitter_1000)[i]))
}
length(indice_ordonne_1000)


donnee_twitter_reorder=donnee_twitter[indice_ordonne,]
donnee_twitter_reorder_1000=donnee_twitter_1000[indice_ordonne_1000,]
color_list= c("#F781BF","#E41A1C","#999999","#377EB8", "#A65628", "#984EA3", "#FF7F00" ,"#FFFF33","#4DAF4A")

##### visualisation des données initiales #####

summary(fblog) # chargement des données
set.seed(42) ### rendre les résulats reproductibles

traitement_donnes=function(verbose){
  fblog2<<-upgrade_graph(fblog)
  party.names<<-sort(unique(V(fblog2)$PolParty)) #### traitement pour future coloration
  party.nums.f<<-as.factor(V(fblog2)$PolParty) ### différents partis
  party.nums<<-as.numeric(party.nums.f) ### numero pour coloration
  fblog2_regroupe <<- contract.vertices(fblog2,party.nums) ### regroupement selon le numero du parti
  E(fblog2_regroupe)$weight <<- 1 #### mettre les poids des edges à 1
  fblog2_simplifie <<- simplify(fblog2_regroupe) ### on enlève les structures complexes
  party.size<<-as.vector(table(V(fblog2)$PolParty)) ## taille du noeud en fonction du nombre de noeuds à l'intérieur 
  if (verbose==T){
    print(table(party.nums.f))  #### distribution des partis parmi les noeuds 
    print(table(V(fblog2)$PolParty)) ### parti politique associé à chaque noeuds 
  }
}


traitement_donnes(verbose=T)


trace_graphe_initiaux=function(type){
  l=layout.kamada.kawai(fblog2) ### pour faire une visualisation claire 
  if (type=="classique"){
    plot(fblog2,layout=l,vertex.label=NA,vertex.color=party.nums,vertex.size=3)  ### premier graph intéressant
  }
  if (type=="simplifie"){
    plot(fblog2_simplifie, vertex.size=5*sqrt(party.size),vertex.label=party.names,vertex.color=color_list,edge.width=sqrt(E(fblog2_simplifie)$weight),vertex.label.dist=1, edge.arrow.size=0)
    }
  else if (type=="total") {
    plot(fblog2_regroupe, vertex.size=5*sqrt(party.size),vertex.label=party.names,vertex.color=color_list,edge.width=sqrt(E(fblog2_regroupe)$weight),vertex.label.dist=0, edge.arrow.size=0)
  }
}
trace_graphe_initiaux("classique")
trace_graphe_initiaux("simplifie")
trace_graphe_initiaux("total")


#### Modèle SBM  
SBM_blog <- mixer(as.matrix(get.adjacency(fblog2)),qmin=12, qmax=12)
model_SBM_blog <- mixer::getModel(SBM_blog) ## permet de récupérer le meilleur modèle \\\ 12 classes ici
noms_variables=SBM_blog$nnames[[1]]
plot(SBM_blog, classes=as.factor(V(fblog2)$PolParty),classes.col=color_list,frame=4)
plot(SBM_blog, classes=as.factor(V(fblog2)$PolParty),frame=4)


prob_predite=SBM_blog$output[[1]]$Taus
prob_predite=as.data.frame(prob_predite)
colnames(prob_predite)<-(1:length(colnames(prob_predite)))





###### les fonctions qui permettent d'avoir les classes des différentes personnalités
classe_appartenance_SBM = function(nb_classe,contenu_classe=FALSE){
  modele = mixer(as.matrix(get.adjacency(fblog2)),qmin=nb_classe,qmax=nb_classe)
  noms_blogs=modele$nnames[[1]]
  list_res=list()
  res_classe = matrix(0,length(colnames(as.data.frame(modele$output[[1]]$Taus))),1)
  for (j in 1:length(colnames(as.data.frame(modele$output[[1]]$Taus)))){
    res_classe[j] <- which.max(modele$output[[1]]$Taus[,j])
  }
  if (contenu_classe == FALSE)
    {
    return(res_classe)
  }
  else{
    for (i in 1:nb_classe){
      list_res[[i]]=noms_blogs[which(res_classe==i)]
      }
    list_res[[nb_classe+1]]=res_classe
    print(table(res_classe))
    return(list_res)
    }
}



classe_appartenance_MBB <- function(graph=as.matrix(get.adjacency(fblog2)), n=8){
  result <- mmsb.collapsed.gibbs.sampler(graph,
                                         K = n, num.iterations=200,
                                         alpha = 0.01,
                                         burnin = 20L,
                                         beta.prior = list(1, diag(5, n) + 1))
  pi <- with(result, t(document_sums) / colSums(document_sums))
  class <- matrix(0,length(V(fblog2)),1)
  for (k in 1:length(V(fblog2))){
    class[k,1]=which.max(pi[k,])
  }
  return(class)
}
classe_appartenance_MBB(as.matrix(get.adjacency(fblog2)),5)

ARI=function(nb_classe,renvoi="indice de Rand"){
  res_classe=classe_appartenance_SBM(nb_classe)
  class=classe_appartenance_MBB(as.matrix(get.adjacency(fblog2)),nb_classe)
  definitif=cbind(res_classe,class)
  colnames(definitif)=c("SBM",'MMB')
  definitif=as.data.frame(definitif)
  tab=table(definitif)
  if (renvoi=="indice de Rand"){
    return(randIndex(tab))
  }
  else if (renvoi=="table"){
    return(tab)}
}

ARI(5)

#### TEST SUR LES ARI #############
comparaison_ARI=function(class_ini,class_fin,renvoi){
  sortie <- lapply(c(class_ini:class_fin), ARI)
  plot(class_ini:class_fin, sortie, "l")
}


comparaison_ARI(5,7)


coord_polygone_regulier=function(n){
  list_x=c()
  list_y=c()
  res=list()
  for (k in (1:n)){
  list_x=c(list_x,cos(2*k*pi/n))
  list_y=c(list_y,sin(2*k*pi/n))
  }
  res[[1]]=list_x
  res[[2]]=list_y
  return(res)
  }

coord_polygone_regulier(5)

vecteur_coordonne=function(n){
  res=coord_polygone_regulier(n)
  matrix_res=matrix(0,n,2)
  for (i in 1:n){
    
    matrix_res[i,]=as.matrix(c(res[[1]][i],res[[2]][i]))
  } 
  return(matrix_res)
}



plot_sur_polygone_SBM=function(n){
  modele<- mixer(as.matrix(get.adjacency(fblog2)),qmin=n, qmax=n)
  matrice_proba=as.matrix(modele$output[[1]]$Taus)
  print(matrice_proba)
  coord=vecteur_coordonne(n)
  res_plot=coord_polygone_regulier(n)
  res_plot[[1]]=c(res_plot[[1]],res_plot[[1]][1])
  res_plot[[2]]=c(res_plot[[2]],res_plot[[2]][1])
  plot(res_plot[[1]],res_plot[[2]],"l",xlab="",ylab="")
  atracer=t(matrice_proba)%*%coord
  vecteur_color=color_list[party.nums]
  points(atracer,col=color_list[party.nums])

}

plot_sur_polygone_SBM(5)


plot_sur_polygone_MMB=function(n){
  result <- mmsb.collapsed.gibbs.sampler (as.matrix(get.adjacency(fblog2)), K=n, num.iterations=100, alpha=0.001, burnin=20L, beta.prior = list( diag(8,n)+1, diag(5,n)+1 ))
  matrice_proba= as.matrix(with(result, t(document_sums) / colSums(document_sums)))
  coord=vecteur_coordonne(n)
  res_plot=coord_polygone_regulier(n)
  res_plot[[1]]=c(res_plot[[1]],res_plot[[1]][1])
  res_plot[[2]]=c(res_plot[[2]],res_plot[[2]][1])
  plot(res_plot[[1]],res_plot[[2]],"l",xlab="",ylab="")
  atracer=matrice_proba%*%coord
  points(atracer,col=color_list[party.nums])
  return(atracer)
}

plot_sur_polygone_MMB(5)


trace_clustering=function(n){
  trace=plot_sur_polygone_MMB(n)
  cl1 = kcca(trace, k=n, kccaFamily("kmeans"))
  image(cl1)
  points(trace, col=color_list[party.nums], pch=19, cex=0.9)
  coord=vecteur_coordonne(n)
  res_plot=coord_polygone_regulier(n)
  res_plot[[1]]=c(res_plot[[1]],res_plot[[1]][1])
  res_plot[[2]]=c(res_plot[[2]],res_plot[[2]][1])
  lines(res_plot[[1]],res_plot[[2]],"l",xlab="",ylab="")
}

trace_clustering(5)




modele=mixer(as.matrix(adja_donnee_twitter),qmin=2,qmax=15)
plot(modele,classes=as.factor(donnee_twitter_reorder$party))


##### par rapport aux données sur Twitter
twitter_color=c("olivedrab","lightblue3","magenta","red2","grey27","navy","royalblue","snow")
twitter_color_1000=c("olivedrab","green","lightblue3","magenta","red2","grey27","navy","royalblue","black","snow")
  

graph_twitter=graph.adjacency(as.matrix(adja_donnee_twitter), mode = "directed")
graph_twitter_1000=graph.adjacency(as.matrix(adja_donnee_twitter_1000), mode = "directed")
graph_twitter=set_vertex_attr(graph_twitter,name="affi_pol",value=as.character(donnee_twitter_reorder$party))
table(donnee_twitter_1000$party)
graph_twitter_1000=set_vertex_attr(graph_twitter_1000,name="affi_pol",value=as.character(donnee_twitter_reorder_1000$party))
affiliation.names<-sort(unique(V(graph_twitter)$affi_pol)) #### traitement pour future coloration
affiliation_1000.names<-sort(unique(V(graph_twitter_1000)$affi_pol)) #### traitement pour future coloration
affiliation.names
affiliation_1000.names
affiliation.nums.f<-as.factor(V(graph_twitter)$affi_pol) ### différents partis
affiliation_1000.nums.f<-as.factor(V(graph_twitter_1000)$affi_pol) ### différents partis
affiliation.nums<-as.numeric(affiliation.nums.f) ### numero pour coloration
affiliation_1000.nums<-as.numeric(affiliation_1000.nums.f) ### numero pour coloration
affiliation.size<-as.vector(table(V(graph_twitter)$affi_pol)) ## taille du noeud en fonction du nombre de noeuds à l'intérieur 
affiliation_1000.size<-as.vector(table(V(graph_twitter_1000)$affi_pol)) ## taille du noeud en fonction du nombre de noeuds à l'intérieur 
affiliation_1000.size
table(affiliation_1000.nums.f)
graph_twitter_regroupe <- contract.vertices(graph_twitter,affiliation.nums) ### regroupement selon le numero du parti
graph_twitter_regroupe_1000 <- contract.vertices(graph_twitter_1000,affiliation_1000.nums)
E(graph_twitter_regroupe)$weight <- 1 #### mettre les poids des edges à 1
E(graph_twitter_regroupe_1000)$weight <- 1 #### mettre les poids des edges à 1
graph_twitter_simplifie <- simplify(graph_twitter_regroupe) ### on enlève les structures complexes
graph_twitter_1000_simplifie <- simplify(graph_twitter_regroupe_1000)

l=layout.kamada.kawai(graph_twitter) ### pour faire une visualisation claire 
l_1000=layout.kamada.kawai(graph_twitter_1000) ### pour faire une visualisation claire 
plot(graph_twitter,layout=l,vertex.label=NA,vertex.color=twitter_color[affiliation.nums],vertex.size=3,edge.arrow.size=.5)
legend(1,1,fill =twitter_color,
       legend=affiliation.names)


plot(graph_twitter_1000,layout=l_1000,vertex.label=NA,vertex.color=twitter_color_1000[affiliation_1000.nums],vertex.size=3,edge.arrow.size=.1)
legend("right",fill =twitter_color_1000,
       legend=affiliation_1000.names,cex=0.5,box.lty=0)
plot(graph_twitter_regroupe, vertex.size=5*sqrt(affiliation.size),vertex.label=affiliation.names,vertex.color=twitter_color,edge.width=sqrt(E(graph_twitter_regroupe)$weight),vertex.label.dist=1, edge.arrow.size=0)
plot(graph_twitter_regroupe_1000, vertex.size=5*sqrt(affiliation_1000.size),vertex.label=affiliation_1000.names,vertex.color=twitter_color_1000,edge.width=sqrt(E(graph_twitter_regroupe_1000)$weight),vertex.label.dist=1, edge.arrow.size=0)
plot(graph_twitter_simplifie, vertex.size=5*sqrt(affiliation.size),vertex.label=affiliation.names,vertex.color=twitter_color,edge.width=sqrt(E(graph_twitter_simplifie)$weight),vertex.label.dist=0, edge.arrow.size=0)
plot(graph_twitter_1000_simplifie, vertex.size=5*sqrt(affiliation_1000.size),vertex.label=affiliation_1000.names,vertex.color=twitter_color_1000,edge.width=sqrt(E(graph_twitter_1000_simplifie)$weight),vertex.label.dist=0, edge.arrow.size=0)































############
data(example)
hergm(d ~ edges_i)

dm <- matrix(1:4, ncol = 2)
dm
set.seed(221)
a=sample(ncol(dm))
a
sample(nrow(dm))
?sample
dm[a,a]



