## importation of the main libraries
library(mixer)
library(ergm)
library(network)
library(igraph)
library("RColorBrewer") ### pour colorer nos graphes
library(flexclust)
library(lda)
library(sand)

setwd(dir="/Users/yoanrussac/desktop")
dir()
buiped=read.csv2(file='biped.csv',sep=",",dec=".")
buiped=as.data.frame(buiped)


color_list= c("#F781BF","#E41A1C","#999999","#377EB8", "#A65628", "#984EA3", "#FF7F00" ,"#FFFF33","#4DAF4A")

##### visualisation des données initiales #####

summary(fblog) # chargement des données
set.seed(42) ### rendre les résulats reproductibles

traitement_donnes=function(verbose){
  fblog=upgrade_graph(fblog)  ### mise à niveau du graphe
  party.names<-sort(unique(V(fblog)$PolParty)) #### traitement pour future coloration
  party.nums.f<-as.factor(V(fblog)$PolParty) ### différents partis
  party.nums=as.numeric(party.nums.f) ### numero pour coloration
  fblog_regroupe <- contract.vertices(fblog,party.nums) ### regroupement selon le numero du parti
  E(fblog_regroupe)$weight <- 1 #### mettre les poids des edges à 1
  fblog_simplifie <- simplify(fblog_regroupe) ### on enlève les structures complexes
  party.size<-as.vector(table(V(fblog)$PolParty)) ## taille du noeud en fonction du nombre de noeuds à l'intérieur 
  if (verbose==T){
    print(table(party.nums.f))  #### distribution des partis parmi les noeuds 
    print(table(V(fblog)$PolParty)) ### parti politique associé à chaque noeuds 
  }
}


trace_graphe_initiaux=function(type){
  l=layout.kamada.kawai(fblog) ### pour faire une visualisation claire 
  if (type=="classique"){
    plot(fblog,layout=l,vertex.label=NA,vertex.color=party.nums,vertex.size=3)  ### premier graph intéressant
  }
  if (type=="simplifie"){
    plot(fblog_simplifie, vertex.size=5*sqrt(party.size),vertex.label=party.names,vertex.color=color_list,edge.width=sqrt(E(fblog_simplifie)$weight),vertex.label.dist=1, edge.arrow.size=0)
    }
  else if (type=="total") {
    plot(fblog_regroupe, vertex.size=5*sqrt(party.size),vertex.label=party.names,vertex.color=color_list,edge.width=sqrt(E(fblog_regroupe)$weight),vertex.label.dist=0, edge.arrow.size=0)
  }
  
}


#### Modèle SBM
SBM_blog <- mixer(as.matrix(get.adjacency(fblog)),qmin=12, qmax=12)
model_SBM_blog <- getModel(SBM_blog) ## permet de récupérer le meilleur modèle \\\ 12 classes ici
plot(SBM_blog, classes=as.factor(V(fblog)$PolParty),classes.col=color_list,frame=4)

prob_predite=fblog.sbm.output$Taus
prob_predite=as.data.frame(prob_predite)
colnames(prob_predite)<-(1:length(colnames(prob_predite)))

noms_variables=fblog.sbm$nnames[[1]]

indice_pue_la_merde=which(buiped$X  %in%  fblog.sbm$nnames[[1]]==FALSE)
buiped=buiped[-indice_pue_la_merde,]
res=matrix(0,length(colnames(prob_predite)),1)
for (j in 1:length(colnames(prob_predite))){
  res[j]<- which.max(prob_predite[,j])
}

table(res)
colnames(res)<- "SBM"
rownames(res)<-noms_variables
res

resultat=cbind(res,buiped)
resultat=resultat[,-2]
table(resultat)





data(blog)

modele=mixer(blog$links,qmin=12,qmax=12)
modele


res_2=matrix(0,length(colnames(modele$output[[1]]$Taus)),1)
for (j in 1:196){
  res_2[j]<- which.max(modele$output[[1]]$Taus[,j])
}
res_2
table(res_2)
table(res)
test_196=cbind(res_2,buiped)
str(test_196)
test_196=test_196[,-2]
colnames(test_196)=c("SBM",'MMB')
tab=table(test_196)
tab
randIndex(tab)
sum(diag(tab))/sum(tab)


classe_appartenance_SBM=function(nb_classe){
  modele=mixer(blog$links,qmin=nb_classe,qmax=nb_classe)
  res_classe=matrix(0,length(colnames(modele$output[[1]]$Taus)),1)
  for (j in 1:196){
    res_classe[j]<- which.max(modele$output[[1]]$Taus[,j])
  }
  return(res_classe)
}

classe_appartenance_MBB <- function(graph, K){
  result <- mmsb.collapsed.gibbs.sampler(graph,
                                         K = K, num.iterations=200,
                                         alpha = 0.01,
                                         burnin = 20L,
                                         beta.prior = list(1, diag(5, K) + 1))
  pi <- with(result, t(document_sums) / colSums(document_sums))
  
  class <- matrix(0,196,1)
  for (k in 1:196){
    class[k,1]=which.max(pi[k,])
  }
  return(class)
  
}


a=MMB(as.matrix(blog$links),4)
table(a)
ARI=function(nb_classe){
  res_classe=classe_appartenance_SBM(nb_classe)
  class=classe_appartenance_MBB(as.matrix(blog$links),nb_classe)
  definitif=cbind(res_classe,class)
  colnames(definitif)=c("SBM",'MMB')
  definitif=as.data.frame(definitif)
  tab=table(definitif)
  return(tab)
  
}

#### TEST SUR LES ARI #############
yu=ARI(11)
yu
line_2=yu[2,]
line_4=yu[4,]
yu2=yu
yu2[2,]<-line_4
yu2[4,]<-line_2


yu_petit=yu[1:3,1:3]
yu_petit=as.table(yu_petit)
yu2_petit=yu_petit
yu_petit
yu2_petit
yu2_petit[1,]<-yu_petit[3,]
yu2_petit[2,]<-yu_petit[1,]
yu2_petit[3,]<-yu_petit[2,]
yu3_petit=yu_petit
yu3_petit[,1]<-yu_petit[,3]
yu3_petit[,2]<-yu_petit[,1]
yu3_petit[,3]<-yu_petit[,2]
yu4_petit=yu3_petit
yu4_petit[1,]<-yu3_petit[3,]
yu4_petit[2,]<-yu3_petit[1,]
yu4_petit[3,]<-yu3_petit[2,]


yu_petit
yu2_petit
yu3_petit
yu4_petit

randIndex(yu_petit)
randIndex(yu2_petit)
randIndex(yu3_petit)
randIndex(yu4_petit)


str(yu)
str(yu_petit)
randIndex(yu)==randIndex(yu2)
sortie <- lapply(c(8:14), ARI)
plot(8:14, sortie, "l")
for (i in 8:13){
  print(i)
  print(ARI(i))
  print("--------")
}
##########################################



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


#######
data(blog)
unique(blog$politicalParty)

fblog



plot_polygone_regulier=function(n){
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

vecteur_coordonne=function(n){
  res=plot_polygone_regulier(n)
  matrix_res=matrix(0,n,2)
  for (i in 1:n){
    
    matrix_res[i,]=as.matrix(c(res[[1]][i],res[[2]][i]))
  } 
  return(matrix_res)
}


plot_sur_polygone_SBM=function(n){
  modele<- mixer(as.matrix(get.adjacency(fblog)),qmin=n, qmax=n)
  matrice_proba=as.matrix(modele$output[[1]]$Taus)
  print(matrice_proba)
  coord=vecteur_coordonne(n)
  res_plot=plot_polygone_regulier(n)
  res_plot[[1]]=c(res_plot[[1]],res_plot[[1]][1])
  res_plot[[2]]=c(res_plot[[2]],res_plot[[2]][1])
  plot(res_plot[[1]],res_plot[[2]],"l",xlab="",ylab="")
  atracer=t(matrice_proba)%*%coord
  vecteur_color=color_list[party.nums]
  points(atracer,col=color_list[party.nums])

}




plot_sur_polygone_MMB=function(n){
  result <- mmsb.collapsed.gibbs.sampler (as.matrix(get.adjacency(fblog)), K=n, num.iterations=100, alpha=0.001, burnin=20L, beta.prior = list( diag(8,K)+1, diag(5,K)+1 ))
  matrice_proba= as.matrix(with(result, t(document_sums) / colSums(document_sums)))
  coord=vecteur_coordonne(n)
  res_plot=plot_polygone_regulier(n)
  res_plot[[1]]=c(res_plot[[1]],res_plot[[1]][1])
  res_plot[[2]]=c(res_plot[[2]],res_plot[[2]][1])
  plot(res_plot[[1]],res_plot[[2]],"l",xlab="",ylab="")
  atracer=matrice_proba%*%coord
  points(atracer,col=color_list[party.nums])
  return(atracer)
}


trace_clustering=function(n){
  trace=plot_sur_polygone_MMB(n)
  cl1 = kcca(trace, k=n, kccaFamily("kmeans"))
  image(cl1)
  points(trace, col=color_list[party.nums], pch=19, cex=0.9)
  coord=vecteur_coordonne(n)
  res_plot=plot_polygone_regulier(n)
  res_plot[[1]]=c(res_plot[[1]],res_plot[[1]][1])
  res_plot[[2]]=c(res_plot[[2]],res_plot[[2]][1])
  lines(res_plot[[1]],res_plot[[2]],"l",xlab="",ylab="")
}











### faire un algo qui pour chaque classe donne la liste des noms de blogs dedans.

