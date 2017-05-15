

#### Dans ce qui suit, nous générons des données simulées sur lesquelles 
#### nous appliquons le MMB


library("lda")
library('mixer')
require(gtools)
require(plotrix)
require(multinomRob)

#we fiw the seed
set.seed(123)

#we generate the graph according to MMSB method

B=matrix(0,4,4)
B[1,1]=0.95
B[2,2]=0.95
B[3,3]=0.4
B[4,4]=0.4

B[1,2]=0.8
B[1,3]=0.3
B[1,4]=0.1

B[2,1]=0.3
B[2,3]=0.8
B[2,4]=0.2

B[3,2]=0.3
B[3,1]=0.1
B[3,4]=0.4

B[4,2]=0.2
B[4,3]=0.8
B[4,1]=0.1

color2D.matplot(B, show.values=2, main="Matrice des interactions")

#Membership vector


mem <- rdirichlet(200, c(0.1,0.1,0.1,0.1))
Z= matrix(0,200,200)
for (i in 1:200){
  for (j in 1:200){
    Z[i,j]=which(rmultinomial(1, mem[i,])==1)
  }
}


#réalisation du graphe
Y=matrix(0,200,200)
for (i in 1:200){
  for (j in 1:200){
    if (runif(1) < B[Z[i,j],Z[j,i]]){
      Y[i,j]=1
    }
  }
}


#on organise Y

reordo <- matrix(0,nrow(Y),1)
for (k in 1:nrow(Y)){
  reordo[k,1]=which.max(mem[k,])
}

verite=Y[order(reordo), order(reordo)]

##on affiche

par(mfrow=c(1,2))
image(t(apply((1-Y), 2, rev)), col= grey.colors(2, start = 0.1, end = 1), xlab = "noeuds", ylab="noeuds")
image(t(apply((1-verite), 2, rev)), col= grey.colors(2, start = 0.1, end = 1), xlab = "noeuds", ylab="noeuds")

## Choix du nombre de clusters

BIC <- function(K){
  result <- mmsb.collapsed.gibbs.sampler(Y,
                                         K = K, num.iterations=300,
                                         alpha = 0.1,
                                         burnin = 20L,
                                         beta.prior = list(diag(8,K)+1, diag(5, K) + 1))
  mem <- with(result, t(document_sums) / colSums(document_sums))
  B <-  with(result, blocks.pos / (blocks.pos + blocks.neg))
  esp <- mem %*% B %*% t(mem)
  proba=(Y*esp)+(1-esp)*(1-Y)
  return(sum(log(proba)) - (K+2*K^2)*(log(ncol(graph))^2))
}

BIC_plot <- lapply(c(2:8), BIC)
par(mfrow=c(1,1))
plot(c(2:8), BIC_plot, "l", xlab="Nombre de classes", ylab="BIC", lwd=2)
grid()

#on choisit 4 clusters



## on utilise le modèle pour alpha =0.01

K=4

result1 <- mmsb.collapsed.gibbs.sampler(Y,
                                       K = K, num.iterations=2000,
                                       alpha = 0.01,
                                       burnin = 20L,
                                       beta.prior = list(diag(8,K)+1, diag(5, K) + 1))

#on retrouve les estimations de B et de pi
memhat1 <- with(result1, t(result1$document_sums) / colSums(result1$document_sums))
Bhat1 <-  with(result1, result1$blocks.pos / (result1$blocks.pos + result1$blocks.neg))

#on réordonne les clusters pour retrouver la matrice initiale
par(mfrow=c(1,1))
color2D.matplot(Bhat1, show.values=TRUE, main="B estimée")

Bhat1ordo <- Bhat1[c(2,3,1,4),c(2,3,1,4)]  #on réordonne à la main les clusters
color2D.matplot(Bhat1ordo, show.values=2, main="B estimée")


#on affiche pi réorganisé

memhat1ordo = memhat1[,c(2,3,1,4)]
ordre1 <- matrix(0,nrow(Y),1)
for (k in 1:nrow(Y)){
  ordre1[k,1]=which.max(memhat1ordo[k,])
}

esp1 = memhat1ordo[order(ordre1),] %*% Bhat1ordo %*% t(memhat1ordo[order(ordre1),])

image(t(apply((1-esp1), 2, rev)), col= grey.colors(5, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds" )




#####   on recommence pour alpha=0.1    ###### 

K=4

result2 <- mmsb.collapsed.gibbs.sampler(Y,
                                        K = K, num.iterations=2000,
                                        alpha = 0.1,
                                        burnin = 20L,
                                        beta.prior = list(diag(8,K)+1, diag(5, K) + 1))

#on retrouve les estimations de B et de pi
memhat2 <- with(result2, t(result2$document_sums) / colSums(result2$document_sums))
Bhat2 <-  with(result2, result2$blocks.pos / (result2$blocks.pos + result2$blocks.neg))

#on réordonne les clusters pour retrouver la matrice initiale
par(mfrow=c(1,1))
color2D.matplot(Bhat2, show.values=2, main="B estimée")

Bhat2ordo <- Bhat2
color2D.matplot(Bhat2ordo, show.values=2, main="B estimée")


#on affiche pi réorganisé

memhat2ordo = memhat2[,c(2,4,3,1)]
ordre2 <- matrix(0,nrow(Y),1)
for (k in 1:nrow(Y)){
  ordre2[k,1]=which.max(memhat2ordo[k,])
}

esp2 = memhat2ordo[order(ordre2),] %*% Bhat2ordo %*% t(memhat2ordo[order(ordre2),])

image(t(apply((1-esp2), 2, rev)), col= grey.colors(5, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds" )


#####   on termine pour alpha=0.25    ######

K=4

result3 <- mmsb.collapsed.gibbs.sampler(Y,
                                        K = K, num.iterations=2000,
                                        alpha = 0.25,
                                        burnin = 20L,
                                        beta.prior = list(diag(8,K)+1, diag(5, K) + 1))

#on retrouve les estimations de B et de pi
memhat3 <- with(result3, t(result3$document_sums) / colSums(result3$document_sums))
Bhat3 <-  with(result3, result3$blocks.pos / (result3$blocks.pos + result3$blocks.neg))

#on réordonne les clusters pour retrouver la matrice initiale
par(mfrow=c(1,1))
color2D.matplot(Bhat3, show.values=TRUE, main="B estimée")

Bhat3ordo <- Bhat3[c(1,4,3,2),c(1,4,3,2)] #coup de chance
color2D.matplot(Bhat3ordo, show.values=2, main="B estimée")


#on affiche pi réorganisé

memhat3ordo = memhat3[,c(1,4,3,2)]
ordre3 <- matrix(0,nrow(Y),1)
for (k in 1:nrow(Y)){
  ordre3[k,1]=which.max(memhat3ordo[k,])
}

esp3 = memhat3ordo[order(ordre3),] %*% Bhat3ordo %*% t(memhat3ordo[order(ordre3),])

image(t(apply((1-esp3), 2, rev)), col= grey.colors(5, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds" )


###Finalement, on affiche les 4 matrices en même temps

par(mfrow=c(2,2))
color2D.matplot(B, show.values=2, main="Matrice des interactions", xlab="", ylab="")
color2D.matplot(Bhat1ordo, show.values=2, main="B estimée (alpha=0.01)", xlab="", ylab="")
color2D.matplot(Bhat2ordo, show.values=2, main="B estimée (alpha=0.1)", xlab="", ylab="")
color2D.matplot(Bhat3ordo, show.values=2, main="B estimée (alpha=0.25)", xlab="", ylab="")

par(mfrow=c(2,2))
image(t(apply((1-verite), 2, rev)), col= grey.colors(2, start = 0.1, end = 1), xlab = "noeuds", ylab="noeuds", main="Graphe")
image(t(apply((1-esp1), 2, rev)), col= grey.colors(5, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds", main="alpha=0.01" )
image(t(apply((1-esp2), 2, rev)), col= grey.colors(5, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds", main="alpha=0.1" )
image(t(apply((1-esp3), 2, rev)), col= grey.colors(5, start = 0.1, end = 1),  xlab = "noeuds", ylab="noeuds", main="alpha=0.25" )



