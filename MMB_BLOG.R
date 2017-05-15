#chargement libraries

library("lda")
library('mixer')
require(gtools)
require(plotrix)

library("ergm")
library("network")
library("igraph")
library("RColorBrewer") ### pour colorer nos graphes
library("flexclust")
library("sand")


##### Premier traitement données ####

summary(fblog)
color_list = c( "#984EA3","#A65628","#4DAF4A","#999999","#FFFF33","#F781BF","#E41A1C","#FF7F00","#377EB8")

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


traitement_donnes(verbose=F)

#### Chargement des données

graph=as.matrix(get.adjacency(fblog2))


##Critère BIC

BIC <- function(K){
  result <- mmsb.collapsed.gibbs.sampler(graph,
                                         K = K, num.iterations=200,
                                         alpha = 0.25,
                                         burnin = 20L,
                                         beta.prior = list(diag(8,K)+1, diag(5, K) + 1))
  mem <- with(result, t(document_sums) / colSums(document_sums))
  B <-  with(result, blocks.pos / (blocks.pos + blocks.neg))
  esp <- mem %*% B %*% t(mem)
  proba=(graph*esp)+(1-esp)*(1-graph)
  return(sum(log(proba)) - (K+2*K^2)*(log(ncol(graph))^2))
}

BIC_simul <- lapply(c(4:15), BIC)
par(mfrow=c(1,1))
plot(c(4:15), BIC_simul, "l", xlab="Nombre de classes", ylab="BIC", lwd=2)
grid()

K=9
result <- mmsb.collapsed.gibbs.sampler(graph,
                                       K = K, num.iterations=1500,
                                       alpha = 0.25, #alpha à faire varier selon les besoins
                                       burnin = 20L,
                                       beta.prior = list(diag(8,K)+1, diag(5, K) + 1))
memb <- with(result, t(document_sums) / colSums(document_sums))
B <-  with(result, blocks.pos / (blocks.pos + blocks.neg))

esp_blog = memb %*% B %*% t(memb)

##on passe à l'affichage


par(mfrow=c(1,3))
image(1-graph, col= grey.colors(2, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds" )
image(1-esp_blog, col= grey.colors(5, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds" )
color2D.matplot(B, show.values=2, xlab="clusters", ylab="clusters")



class <- matrix(0,ncol(graph),1)  #on récupère la classe de chacun des blogs
for (k in 1:ncol(graph)){
  class[k,1]=which.max(memb[k,])
}

B_corrige=B  #matrice corrigée
for (i in 1:K){
  for (j in 1:K){
    B_corrige[i,j]=B[i,j]*sum(class==i)*sum(class==j)
  }
}

par(mfrow=c(1,1))
color2D.matplot(B_corrige, show.values=1, xlab="clusters", ylab="clusters")

list_partis=list() #composition des classes en termes de partis
for (i in 1:K){
  list_partis[[i]]=party.names[party.nums[which(class==i)]]
}

par(mfrow=c(3,3)) #on affiche les histogrammes des compositions
for (k in 1:K){
  plot(table(list_partis[k]), main=paste("Classe numéro ",k), xlab="partis", ylab="nombre de blogs", lwd=6 )
}

## on projette le modèle estimé dans un polygone

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


vecteur_coordonne=function(n){
  res=coord_polygone_regulier(n)
  matrix_res=matrix(0,n,2)
  for (i in 1:n){
    
    matrix_res[i,]=as.matrix(c(res[[1]][i],res[[2]][i]))
  } 
  return(matrix_res)
}

plot_sur_polygone=function(n, matrice_proba){
  coord=vecteur_coordonne(n)
  res_plot=coord_polygone_regulier(n)
  res_plot[[1]]=c(res_plot[[1]],res_plot[[1]][1])
  res_plot[[2]]=c(res_plot[[2]],res_plot[[2]][1])
  vecteur_color=color_list[party.nums]
  coord=vecteur_coordonne(n)
  atracer=matrice_proba%*%coord
  
  cl1 = kcca(atracer, k=n, kccaFamily("kmeans"))  #on utilise les k-means pour créer des classes
                                                  # plus pertinentes
  image(cl1)
  par(new=TRUE)
  lines(res_plot[[1]],res_plot[[2]],"l",xlab="",ylab="", lwd=3, xlim=c(-1.3, 1.1), ylim=c(-1.1, 1.3) )
  for (m in 1:n){
    lines(c(0,res_plot[[1]][m]), c(0,res_plot[[2]][m]), "l", lty=2, lwd=0.9, col="gray90")
  }
  points(atracer,col=color_list[party.nums], pch=20)
  
  centroids=matrix(0,length(list_partis),n)   # affiche le barycentre de chaque parti
  for (k in 1:length(party.names)){
    centroids[k,]=colMeans(matrice_proba[ party.nums == k , ])
  }
  coord_centro= centroids %*% coord
  points(coord_centro, col=color_list[], pch=15, cex=2)
  
}

par(mfrow=c(1,1))
plot_sur_polygone(9, memb)
legend(x=-1.35, y=1.25, legend = party.names, col=color_list, pch=20, text.col = color_list, text.width = 0.25, x.intersp	
=0.5, cex=0.7, bg="gray90")



### On calcule l'ARI

definitif=cbind(class,party.nums)
colnames(definitif)=c("MMB",'Parti')
definitif=as.data.frame(definitif)
tab=table(definitif)
randIndex(tab) #nous renvoie l'ARI
#32.36%
### on généralise la procédure précédente 

classe_appartenance_MBB <- function(graph=as.matrix(get.adjacency(fblog2)), n=8){
  result <- mmsb.collapsed.gibbs.sampler(graph,
                                         K = n, num.iterations=200,
                                         alpha = 0.1,
                                         burnin = 20L,
                                         beta.prior = list(1, diag(5, n) + 1))
  pi <- with(result, t(document_sums) / colSums(document_sums))
  class <- matrix(0,length(V(fblog2)),1)
  for (k in 1:length(V(fblog2))){
    class[k,1]=which.max(pi[k,])
  }
  return(class)
}

ARI=function(nb_classe,renvoi="indice de Rand"){
  res_classe=party.nums
  class=classe_appartenance_MBB(as.matrix(get.adjacency(fblog2)),nb_classe)
  definitif=cbind(res_classe,class)
  colnames(definitif)=c("Partis",'MMB')
  definitif=as.data.frame(definitif)
  tab=table(definitif)
  if (renvoi=="indice de Rand"){
    return(randIndex(tab))
  }
  else if (renvoi=="table"){
    return(tab)}
}


comparaison_ARI_partis=function(class_ini,class_fin,renvoi){
  sortie <- lapply(c(class_ini:class_fin), ARI)
  plot(class_ini:class_fin, sortie, "l", xlab="nombre de clusters", ylab="ARI", lwd=2)
  grid()
}


comparaison_ARI_partis(5,15)  #permet d'afficher l'ARI en fonction du nombre de classes



#### Annexe: convergence de la log-vraisemblance en fonction du nombre d'itérations

conv <- function(n){
  result <- mmsb.collapsed.gibbs.sampler(graph,
                                         K = K, num.iterations=n,
                                         alpha = 0.25,
                                         burnin = 20L,
                                         beta.prior = list(diag(8,K)+1, diag(5, K) + 1))
  mem <- with(result, t(document_sums) / colSums(document_sums))
  B <-  with(result, blocks.pos / (blocks.pos + blocks.neg))
  esp <- mem %*% B %*% t(mem)
  proba=(graph*esp)+(1-esp)*(1-graph)
  return(sum(log(proba)))
}

test=lapply(50*c(1:20), conv)

par(mfrow=c(1,1))
plot(50*c(1:10), test, "l", xlab="nombre d'itérations", ylab= "log-vraisemblance", lwd=2)
grid()


