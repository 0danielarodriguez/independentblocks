library(mvtnorm)
library(ggplot2)
library(mclust)

#This function performs cross-validation (CV) and penalty selection using cross-validation, with the same folds
# data in  a nxp matrix 
# divergence to compute the estimator
# distance to compute cross validation. 
# lambda_max: maximum value for choosing lambda
# step: for constructing a grid for lambda
# k: for k-folder

eip_final <- function(data,lambda_max=2, step=0.01,k=5,divergencia=divergencia_Frobenius, distancia=distancia_Frobenius,correlation=FALSE){
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  shuffle <- sample(1:n,n,replace = FALSE)
  grilla_lambda <- seq(0,lambda_max,by=step)
  loss_lambda <- rep(NA, length(grilla_lambda))
  
  bs= floor(n/k) # for cross validation
  
  for(i in 1:length(grilla_lambda)){
    lambda <- grilla_lambda[i]
    cum_loss <- 0
    for (j in 1:k)
    {
      test_index <- shuffle[seq((j-1)*bs+1,j*bs)]
      train_index <-shuffle[-seq((j-1)*bs+1,j*bs)]
      
      X_test <- data[test_index,]
      X_train <- data[train_index,]  
      
      eip_with_lambda(X_train,lambda, divergencia)$Sigma_est_part
      
      if(correlation==TRUE){
        l_lambda <-distancia( matriz_correlacion(eip_with_lambda(X_train,lambda, divergencia)$Sigma_est_part),matriz_correlacion(cov(X_test)))
      }
      
      if(correlation==FALSE){
        l_lambda <-distancia(eip_with_lambda(X_train,lambda, divergencia)$Sigma_est_part,cov(X_test))
      }
    
      cum_loss <- cum_loss+l_lambda
    }
    loss_lambda[i] <- cum_loss
  }
  
  lambda_opt <- grilla_lambda[which.min(loss_lambda)]
  
  
  Sigma_est <- cov(data)
  list_of_partitions <- all_partitions(seq(1:p))
  
  final_part_id<- eip_with_lambda(data, lambda_opt, divergencia)$est_part_id
  final_partition <-  list_of_partitions[[final_part_id]]
  final_Sigma_est_part <- varianza_en_bloques(Sigma_est,final_partition)
  salida_1 <- list(est_part_id=final_part_id,estimated_partition=final_partition,
                   Sigma_est_part=final_Sigma_est_part, lambda_opt=lambda_opt)      
  
  partitions <- all_partitions(seq(1:p))
  loss_en_partitions <- rep(NA, length(partitions))
  
  bs= floor(n/k) # for cross validation
  
  for(i in 1:length(partitions)){
    partition <- partitions[[i]]
    cum_loss <- 0
    for (j in 1:k)
    {
      test_index <- shuffle[seq((j-1)*bs+1,j*bs)]
      train_index <-shuffle[-seq((j-1)*bs+1,j*bs)]
      
      X_test <- data[test_index,]
      X_train <- data[train_index,]  
      
      Sigma_test <- cov(X_test)
      Sigma_train <- varianza_en_bloques(cov(X_train), partition)
      
      
      if(correlation==TRUE)
      {
        Sigma_test <- matriz_correlacion(Sigma_test)
        Sigma_train <- matriz_correlacion(Sigma_train)
        
      }
      
      
      l_partition <- distancia(Sigma_test,Sigma_train)
      
      cum_loss <- cum_loss+l_partition
    }
    loss_en_partitions[i] <- cum_loss
  }
  
  optimal_partition_id <- which.min(loss_en_partitions)
  optimal_partition <- partitions[[optimal_partition_id]]
  
  
  Sigma_est <- cov(data)
  final_Sigma_est_part <- varianza_en_bloques(Sigma_est,optimal_partition)
  salida_2 <- list(est_part_id_cv=optimal_partition_id,estimated_partition_cv=optimal_partition,
                   Sigma_est_part_cv=final_Sigma_est_part)      
  salida <- c(salida_1,salida_2)
  
  salida
}



# Final Alternative Functions
#eip: Estimated Independent Partition

#data in  a nxp matrix 
#lambda is the penalization factor.
# divergencia for comparing matrices
# correlation==TRUE use correlation matrices for defining the estimator. 
#returns a list with: 
# est_part_id: estimated partition id
# estimated_partition: partition minimizing the penalized loss with lambda
# Sigma_est_part: Estimated covariance matrix under independence with estimated partition

eip_with_lambda <- function(data, lambda, divergencia=divergencia_Frobenius, correlation=FALSE){
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  list_of_partitions <- all_partitions(seq(1:p))
  Sigma_est <- cov(data)
  everything <- matrix(nrow = length(list_of_partitions),ncol=3)# to store part_id, loss and part_size
  for (j in 1:length(list_of_partitions)){
  
    if(correlation==TRUE){
      matriz_estimada <- matriz_correlacion(Sigma_est)
    }
    
    if(correlation==FALSE){
      matriz_estimada <- Sigma_est
    }
    
    everything[j,] <- c(j, divergencia(matriz_estimada, list_of_partitions[[j]]), 
                        length(list_of_partitions[[j]]))
  }
  everything <- as.data.frame(everything)
  colnames(everything) <- c("part_id", "loss", "part_size" )
  
  
  estimated_part_id <-  everything$part_id[which.min(everything$loss+lambda/everything$part_size)]
  estimated_partition <-  list_of_partitions[[estimated_part_id]]
  Sigma_est_part <- varianza_en_bloques(Sigma_est,estimated_partition)
  salida <- list(est_part_id=estimated_part_id,estimated_partition=estimated_partition,Sigma_est_part=Sigma_est_part)      
  salida
  
}

# to the cross validation procedure
#returns the partition minimizing the penalized loss  and lambda_opt
# data in  a nxp matrix 
# divergencia para calcuar el estimador 
# distancia para hacer cross validation. 
# lambda_max: maximum value for choosing lambda
# step: for constructing a grid for lambda
# k: for k-folder

eip_cv <- function(data,lambda_max=2, step=0.01,k=5,divergencia=divergencia_Frobenius, distancia=distancia_Frobenius,correlation=FALSE){
  n <- dim(data)[1]
  p <- dim(data)[2]
  shuffle <- sample(1:n,n,replace = FALSE)
  grilla_lambda <- seq(0,lambda_max,by=step)
  loss_lambda <- rep(NA, length(grilla_lambda))
  
  bs= floor(n/k)
  
  for(i in 1:length(grilla_lambda)){
    lambda <- grilla_lambda[i]
    cum_loss <- 0
    for (j in 1:k)
    {
      test_index <- shuffle[seq((j-1)*bs+1,j*bs)]
      train_index <-shuffle[-seq((j-1)*bs+1,j*bs)]
      
      X_test <- data[test_index,]
      X_train <- data[train_index,]  
      
      eip_with_lambda(X_train,lambda, divergencia)$Sigma_est_part
      
      if(correlation==TRUE){
        l_lambda <-distancia( matriz_correlacion(eip_with_lambda(X_train,lambda, divergencia)$Sigma_est_part),matriz_correlacion(cov(X_test)))
      }
      
      if(correlation==FALSE){
        l_lambda <-distancia(eip_with_lambda(X_train,lambda, divergencia)$Sigma_est_part,cov(X_test))
      }
      
      cum_loss <- cum_loss+l_lambda
    }
    loss_lambda[i] <- cum_loss
  }
  
  lambda_opt <- grilla_lambda[which.min(loss_lambda)]
  
  
  Sigma_est <- cov(data)
  list_of_partitions <- all_partitions(seq(1:p))
  
  final_part_id<- eip_with_lambda(data, lambda_opt, divergencia)$est_part_id
  final_partition <-  list_of_partitions[[final_part_id]]
  final_Sigma_est_part <- varianza_en_bloques(Sigma_est,final_partition)
  salida <- list(est_part_id=final_part_id,estimated_partition=final_partition,
                 Sigma_est_part=final_Sigma_est_part, lambda_opt=lambda_opt)      
  salida
}



# Compute the partition without penalization

cv_traditional <- function(data, step=0.01,k=5, distancia=distancia_Frobenius,correlation=FALSE){
  n <- dim(data)[1]
  p <- dim(data)[2]
  shuffle <- sample(1:n,n,replace = FALSE)
  partitions <- all_partitions(seq(1:p))
  loss_en_partitions <- rep(NA, length(partitions))
  
  bs= floor(n/k) 
  
  for(i in 1:length(partitions)){
    partition <- partitions[[i]]
    cum_loss <- 0
    for (j in 1:k)
    {
      test_index <- shuffle[seq((j-1)*bs+1,j*bs)]
      train_index <-shuffle[-seq((j-1)*bs+1,j*bs)]
      
      X_test <- data[test_index,]
      X_train <- data[train_index,]  
      
      Sigma_test <- cov(X_test)
      Sigma_train <- varianza_en_bloques(cov(X_train), partition)
      
      
      if(correlation==TRUE)
      {
        Sigma_test <- matriz_correlacion(Sigma_test)
        Sigma_train <- matriz_correlacion(Sigma_train)
        
      }
      
      
      l_partition <- distancia(Sigma_test,Sigma_train)
      
      cum_loss <- cum_loss+l_partition
    }
    loss_en_partitions[i] <- cum_loss
  }
  
  optimal_partition_id <- which.min(loss_en_partitions)
  optimal_partition <- partitions[[optimal_partition_id]]
  
  
  Sigma_est <- cov(data)
  final_Sigma_est_part <- varianza_en_bloques(Sigma_est,optimal_partition)
  salida <- list(est_part_id=optimal_partition_id,estimated_partition=optimal_partition,
                 Sigma_est_part=final_Sigma_est_part)      
  salida
}





###############################
#Auxiliar Functions

#return a pxp matrix with 1 on the positions  (i,j), with i,j in  vec, 
#and 0 otherwise
bloque <- function(p,vec)
{
 M <- matrix(0,p,p) 
 for(i in 1:length(vec))
 {
   for (j in 1:length(vec))
     
    M[vec[i],vec[j]] <- 1 
 }
M
}

# compute Sigma from blocks

varianza_en_bloques <- function(Sigma,particion)
{
p <- dim(Sigma)[1]
salida <- matrix(0,p,p) 
for (i in 1:length(particion))
{
 salida <- salida+bloque(p,particion[[i]]) 
}
salida <- salida*Sigma
salida
}


#Matrixs Distnaces

distancia_Frobenius <- function(Sigma_1,Sigma_2) # es el cuadrado
  {
  salida <- sum((Sigma_1-Sigma_2)^2)
  salida
}

distancia_Hellinger <- function(Sigma_1,Sigma_2, ...)
  {
resto <- det(Sigma_1)^{1/4}*det(Sigma_2)^{1/4}/det((Sigma_1+Sigma_2)/2)^{1/2}
medias <- list(...)
if (length(medias)>0 ){
    factor <- exp((-1/8)*(medias[[1]]-medias[[2]])%*%solve((Sigma_1+Sigma_2)/2)%*%
        (medias[[1]]-medias[[2]]))
resto <-resto*factor
}
salida <- 1-resto
salida}



distancia_KL <- function(Sigma_1,Sigma_2,...) 
{  p <- dim(Sigma_1)[1]
  Sigma_1_mu <- solve(Sigma_1)
  salida <-  sum(diag(Sigma_1_mu%*%Sigma_2)) -log(det(Sigma_1_mu%*%Sigma_2))-p

  medias <- list(...)
  if (length(medias)>0 ){
    dif_medias <- (medias[[1]]-medias[[2]])%*%solve(Sigma_1)%*%
      (medias[[1]]-medias[[2]])
    salida <- salida+dif_medias
  }
  
  salida/2
  }


distancia_W <- function(Sigma_1,Sigma_2)
{  descompongo1 <- svd(Sigma_1)
  raiz_Sigma_1 <- descompongo1$u%*%diag((descompongo1$d)^{1/2})%*%t(descompongo1$v)
  
  grande <- raiz_Sigma_1%*%Sigma_2%*%raiz_Sigma_1

  descompongogrande <- svd(grande)
  raiz_grande <- descompongogrande$u%*%diag((descompongogrande$d)^{1/2})%*%t(descompongogrande$v)
  
  salida <-sum(diag(Sigma_1+Sigma_2-2*raiz_grande))
  
    
  salida
}

#Compute the Divergence or disntace between   Sigma and pi-Sigma, 

divergencia_Frobenius <- function(Sigma,particion){
  Sigma_partida <- varianza_en_bloques(Sigma,particion)
  salida <- distancia_Frobenius(Sigma,Sigma_partida)
  salida
}


divergencia_Hellinger <- function(Sigma,particion){
  Sigma_partida <- varianza_en_bloques(Sigma,particion)
  salida <- distancia_Hellinger(Sigma,Sigma_partida)
  salida
}


divergencia_KL <- function(Sigma,particion){
  Sigma_partida <- varianza_en_bloques(Sigma,particion)
  salida <- distancia_KL(Sigma,Sigma_partida)
  salida
}


divergencia_W <- function(Sigma,particion){
  Sigma_partida <- varianza_en_bloques(Sigma,particion)
  salida <- distancia_W(Sigma,Sigma_partida)
  salida
}



matriz_correlacion <- function(Sigma){
  p <- dim(Sigma)[1]
  salida <- matrix(NA,p,p)
  for(j in 1:p){
    for(i in 1:p){
      salida[i,j] <- Sigma[i,j]/sqrt(Sigma[i,i]*Sigma[j,j])
    }
  }
  salida
}

#compute all the partitions

all_partitions <- function(set){
  todas=list(list())
  
  for (i in 1:length(set)){
    elemento=set[i]
    nuevas=list()
    for(j in 1:length(todas)){
      
      nuevas=c(nuevas,  mistura(todas[[j]],elemento))
      
    }
    todas <- nuevas
  }
  
  todas
  
}

mistura <- function(particion, elemento) { 
  if(length(particion) == 0){
    agrego <-c(elemento)  
    salida <- list(agrego)
  }
  if (length(particion) > 0) {

    agrego <- c(particion,c(elemento)) 
    salida <- list(agrego) 
    for (i in 1:length(particion)) {
      
      viejo=particion[-i] 
      pp=c(particion[[i]], elemento) 
      nueva=c(viejo,list(pp))
      salida <- c(salida, list(nueva))
    }
  }
  

  return(salida)
}

