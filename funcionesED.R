#indice de validez
#funcion que calcula el fitness del individuo
homogeneidad <- function(datos, individuo, numCluster){
  #calcular la matriz de distancia
  matrizPearson <- matPearson(datos)
  promHomogeneidad <- NULL
  resultadoCluster <- NULL
  distancia <- NULL
  
  for (i in 1:numCluster) {
    #cantidad de elementos que tiene un cluster
    items <- length(individuo[[2]][[i]])
    
    #si items tiene valor 1 el calculo de homogeneidad es 1
    if(items == 1){
      promHomogeneidad <- c(promHomogeneidad, 1)
    }
    
    else{
      for (j in 1:items) {
        
        #obtener dsitancia promedio del item a los demas items del cluster
        for(k in  1:items){
          if(j != k){
            distancia <- c(distancia, matrizPearson[individuo[[2]][[i]][[j]], individuo[[2]][[i]][[k]]])
          }
        }
        
        resultado <- mean(distancia) / ((items * (items - 1)) / 2)
        resultadoCluster <- c(resultadoCluster, resultado)
        distancia <- NULL
      }
    }
    #se guarda el resultado del cluster
    promHomogeneidad <- c(promHomogeneidad, sum(resultadoCluster))
    resultadoCluster <- NULL
  }
  
  return(sum(promHomogeneidad)/numCluster)
}

#funcion para crear la matriz de distancia
matPearson <- function(datos){
  #se crea matriz de distancia con dimensiones igual al numero de genes que tiene los datos
  dimension <- dim(datos)
  matrizPearson <- diag(1, dimension[2])
  
  #se calcula la matriz de distancia
  for(i in 1:(dimension[2] - 1)){
    for (j in (i + 1):dimension[2]) {
      matrizPearson[i,j] <- 1 - pearson(datos[,i], datos[,j])
      matrizPearson[j,i] <- matrizPearson[i,j]
    }
  }
  
  return(matrizPearson)
}

#Funcion de distancia 
#Correlacion de pearson entre dos vectores
pearson <- function(i, j){
  if(var(j) == 0){
    return(0)
  }
  x <- cov(i,j)/((var(i) * var(j))^(0.5))
  return(x)
}

#Funcion que crea la poblacion inicial 
setup <- function(datos, p, numCluster, features , tamano, cantDatos){
  #tamaño de la poblacion, n entero , features cantidad de caracteristicas
  individuo <-NULL
  #p es la poblacion inicial con largo igual al parametro tamaño
  p <- list()
  p <- vector("list", length = tamano)
  
  #se crean los individuos y se guardan en p
  for(i in 1:tamano){
    individuo <- generarIndividuo(datos, numCluster, features, cantDatos)
    p[[i]] <-  individuo
  }
  return(p)
}

#Funcion que genera los individuos 
generarIndividuo <- function(datos, numCluster, features, cantDatos){
  #se crea la matriz de centroides 
  centroides <- numeric(numCluster * features) 
  dim(centroides) <- c(numCluster, features)
  
  #se eligen cetroides aleatorios iniciales 
  #centroidesRandom <- sample(1:dimensionesDatos[2], numCluster, replace = FALSE)
  centroidesRandom <- sample(3:cantDatos, numCluster, replace = FALSE)
  for (i in 1:numCluster) {
    for (j in 1:features) {
      centroides[i,j] <- datos [j,centroidesRandom[i]]
    }
  }
  #Se crean los cluster como una lista vacia
  clusters <- list()
  #se crea el individuo compuesto por los centroides y clusters
  individuo <- list(centroides, clusters)
  #Se generan clusters del individuo 
  individuo <- ComputeClustering(individuo, datos, numCluster, cantDatos)
  return(individuo)
}

#Funcion que obtiene el centroide mas cercano a un item
calcularDistancia <- function(dato, centroides, numCluster){
  #se el centroide con menor distancia a un item  
  distanciaMenor <- pearson(dato, centroides[1,])
  posCentroide <- 1
  for (i in 2:numCluster) {
    distancia <- pearson(dato, centroides[i,])
    if(distanciaMenor == 1){
      resultado <- c(0, i)
      return(resultado)
    }
    if(distancia > distanciaMenor){
      distanciaMenor <- distancia
      posCentroide <- i
    }
    
  }
  resultado <- c(1, posCentroide)
  return(resultado)
}

#Funcion que calcula el centroide de un cluster
calcularCentroide <- function(datos, individuo, cantidadItem){
  #si el cluster tiene sol un item, se retorna el individuo
  if(cantidadItem == 1){
    return(individuo)
  }
  
  nuevoCentroide <- datos[,individuo[1]]
  for (i in 2:cantidadItem) {
    nuevoCentroide <- nuevoCentroide + datos[,individuo[i]]
  }
  nuevoCentroide <- nuevoCentroide/cantidadItem
  return(nuevoCentroide)
  
}

#Funcion que genera los cluster de los individuos
ComputeClustering <- function(individuo, datos, numCluster, cantDatos){
  
  #lista de cluster vacio
  individuo[[2]] <- vector("list", length = numCluster)
  
  for(i in 1:cantDatos){
    #Se busca el centroide mas cercano al item
    ClusterCerca <- calcularDistancia(datos[,i], individuo[[1]], numCluster)
    
    #Si cluster cerca es 0, es por que es el centroide
    if(ClusterCerca[1] == 0){
      individuo[[2]][[ClusterCerca[2]]] <- c(individuo[[2]][[ClusterCerca[2]]], i)
    }
    else{
      #añadir el item al cluster mas cercano
      individuo[[2]][[ClusterCerca[2]]] <- c(individuo[[2]][[ClusterCerca[2]]], i)
      
    }
  }
  
  nuevoCentroide <- NULL
  
  #Se calcula el nuevo centroide del cluster
  for(j in 1:numCluster){
    cantidadItem <- length(individuo[[2]][[j]])
    if(cantidadItem != 0){
      nuevoCentroide <- calcularCentroide(datos, individuo[[2]][[j]], cantidadItem)
      individuo[[1]][j,] <- nuevoCentroide
    }
  }
  
  return(individuo)
}

#Funcion que mueve los centroides vacios al centroide no vacio mas cercano
buscarCentroides <- function(centroides, posCentroide, clusters, numCluster){
  distanciaMin <- -1
  posMin <- 0
  for (i in 1:numCluster) {
    if(i != posCentroide && length(clusters[[i]] > 1)){
      distancia <- pearson(centroides[posCentroide,], centroides[i,])
      if(distancia > distanciaMin){
        distanciaMin <- distancia
        posMin <- i
      }
    }
  }
  
  return(posMin)
}

#Funcion que realiza las etapas de mutacion y cruce del algoritmo
crossover <- function(p, a, b, c, datos, fMutacion, CR, numCluster, features, cantDatos){
  #Entrada: p individuo padre, w y z individuos elegidos aleatoriamente, g numero de generaciones, d medida de distancia
  #datos dataset, fMutacion es el factor de mutacion, CR es el cruce
  
  #crear vector de cluster de prueba
  clusterPrueba <- list()
  
  #crea la matriz de centroides de prueba
  centroidesPrueba <- numeric(numCluster * features) 
  dim(centroidesPrueba) <- c(numCluster, features)
  
  for (i in 1:numCluster) {
    if(runif(1, 0.0, 1.0) < CR || i== 4 ){
      for (j in 1:features) {
        centroidesPrueba[i,j] <- a[[1]][i,j] + (fMutacion * (b[[1]][i,j] - c[[1]][i,j])) 
      }
    }else{
      centroidesPrueba[i,] <- p[[1]][i,]
    }
  }
  #aca hacer la condicion
  
  #Se crea el individuo hijo con los centroide y el cluster de prueba
  iHijo <- list(centroidesPrueba, clusterPrueba)
  
  #Se generan clusters del individuo hijo
  individuoPrueba <- ComputeClustering(iHijo, datos, numCluster, cantDatos)
  
  contador <- 0
  while (1) {
    
    #Se revisa que ningun cluster este vacio
    for (i in 1:numCluster) {
      if(length(individuoPrueba[[2]][[i]]) == 0){
        posCentroide <- buscarCentroides(individuoPrueba[[1]], i, individuoPrueba[[2]], numCluster)
        centroidesPrueba[i,] <- (centroidesPrueba[i,] + centroidesPrueba[posCentroide,]) / 2
        contador <- contador + 1
      }
    }
    
    if(contador == 0){
      break
    }
    
    #Se generan clusters del individuo hijo
    clusterPrueba <- list()
    iHijo <- list(centroidesPrueba, clusterPrueba)
    individuoPrueba <- ComputeClustering(iHijo, datos, numCluster, cantDatos)
    contador <- 0
    
  }
  
  return(individuoPrueba)
}

#Funcion que compara el fitness del individuo padre y el individuo hijo o de prueba
reemplazo <- function(datos, iPadre, iPrueba, features, numCluster){
  
  homogeneidadPadre <- homogeneidad(datos, iPadre, numCluster) 
  homogeneidadHijo <- homogeneidad(datos, iPrueba, numCluster)
  
  if(homogeneidadPadre < homogeneidadHijo){
    return(iPadre)
  }else{
    return(iPrueba)
  }
}

#Funcion que obtiene el mejor individuo de la poblacion
bestIndividuo <- function(datos, poblacion, features, numCluster, tamano){
  
  mejorIndividuo <- homogeneidad(datos, poblacion[[1]], numCluster)
  posicion <- 1
  for (i in 2:tamano) {
    individuo <- homogeneidad(datos, poblacion[[i]], numCluster)
    
    if(individuo < mejorIndividuo){
      mejorIndividuo <- individuo
      posicion <- i
    }
  }
  result <- c(mejorIndividuo, posicion)
  return(result)
}

#Funcion que obtiene el fitness de todos los individuos de la poblacion
individuosRend <- function(datos, poblacion, numCluster, tamano){
  
  fitness <- NULL
  for (i in 1:tamano) {
    fitness <- c(fitness, homogeneidad(datos, poblacion[[i]], numCluster))
  }
  return(fitness)
}