## code to prepare `DATASET` dataset goes here

##--------------------------------------------------------------------------
## Simulate Data for two-component
##--------------------------------------------------------------------------
simulate_2comp <- function(G = 500, My = 100, M1 = 100,
                           output.more.info = FALSE){
  requireNamespace("truncdist", quietly=TRUE)
  requireNamespace("SummarizedExperiment", quietly=TRUE)
  ## Simulate MuN and MuT for each gene
  MuN <- rnorm(G, 7, 1.5)
  MuT <- rnorm(G, 7, 1.5)
  Mu <- cbind(MuN, MuT)
  colnames(Mu) <- c('MuN', 'MuT')
  rownames(Mu) <- paste('Gene', seq = seq(1, G))
  
  ## Simulate SigmaN and SigmaT for each gene
  SigmaN <- runif(n = G, min = 0.1, max = 0.8)
  SigmaT <- runif(n = G, min = 0.1, max = 0.8)
  Sigma <- cbind(SigmaN, SigmaT)
  colnames(Sigma) <- c('SigmaN', 'SigmaT')
  rownames(Sigma) <- paste('Gene', seq = seq(1, G))
  
  ## Initial values
  data.N1 <- matrix(0, G, M1)
  data.Y <- matrix(0, G, My)
  True.data.N1 <- matrix(0, G, My)
  True.data.T <- matrix(0, G, My)
  ## Creat row and column name
  rownames(data.N1) <- rownames(data.Y) <- 
    rownames(True.data.N1) <- rownames(True.data.T) <-
    paste('Gene', seq = seq(1, G))
  colnames(data.N1) <- paste('Sample', seq = seq(1, M1))
  colnames(data.Y) <- colnames(True.data.N1) <- 
    colnames(True.data.T) <- paste('Sample', seq = seq(1, My))
  
  ## Simulate Tumor Proportion
  PiT = truncdist::rtrunc(n = My,
                          spec = 'norm', 
                          mean = 0.55,
                          sd = 0.2,
                          a = 0.25,
                          b = 0.95)
  pi <- rbind(1-PiT, PiT)
  rownames(pi) <- c('PiN', 'PiT')
  colnames(pi) <- paste('Sample', seq = seq(1, My))
  
  ## Simulate Data
  for(k in 1:G){
    
    data.N1[k,] <- 2^rnorm(M1, MuN[k], SigmaN[k]); # normal reference
    
    True.data.T[k,] <- 2^rnorm(My, MuT[k], SigmaT[k]);  # True Tumor
    
    True.data.N1[k,] <- 2^rnorm(My, MuN[k], SigmaN[k]);  # True Normal
    
    data.Y[k,] <- pi[1,]*True.data.N1[k,] + pi[2,]*True.data.T[k,] # Mixture Tumor
    
  }
  
  ## Transfer into bioconductor format
  data.Y <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(data.Y)),
                                 rowData = DataFrame(row.names = rownames(data.Y)),
                                 colData = DataFrame(row.names = colnames(data.Y)))
  data.N1 <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(data.N1)),
                                  rowData = DataFrame(row.names = rownames(data.N1)),
                                  colData = DataFrame(row.names = colnames(data.N1)))
  True.data.T <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(True.data.T)),
                                      rowData = DataFrame(row.names = rownames(True.data.T)),
                                      colData = DataFrame(row.names = colnames(True.data.T)))
  True.data.N1 <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(True.data.N1)),
                                       rowData = DataFrame(row.names = rownames(True.data.N1)),
                                       colData = DataFrame(row.names = colnames(True.data.N1)))
  
  if(output.more.info){
    test.data = list(pi=pi, Mu = Mu, Sigma = Sigma,
                     data.Y = data.Y, data.N1 = data.N1,
                     True.data.T = True.data.T, True.data.N1 = True.data.N1)
  }else{
    test.data = list(pi=pi, Mu = Mu, Sigma = Sigma,
                     data.Y = data.Y, data.N1 = data.N1)
  }
  
  return(test.data)
}
set.seed(123)
test.data.2comp = simulate_2comp(G = 500, My = 100, M1 = 100)
usethis::use_data(test.data.2comp, compress = 'xz', overwrite = TRUE)



##--------------------------------------------------------------------------
## Simulate Data for three-component mixed cell line
##--------------------------------------------------------------------------
simulate_3comp <- function(G1 = 475, G2 = 25, My = 5, M1 = 100, M2 = 100,
                           output.more.info = FALSE){
  requireNamespace("truncdist", quietly=TRUE)
  requireNamespace("SummarizedExperiment", quietly=TRUE)
  G = G1 + G2
  ## Simulate MuN1, MuN2 and MuT for each gene
  MuN1 <- rnorm(G1 + G2, 7, 1.5)
  MuN2_1st <- MuN1[1:G1] + truncdist::rtrunc(n = 1, 
                                             spec = 'norm',
                                             mean = 0,
                                             sd = 1.5,
                                             a = -0.1, 
                                             b = 0.1)
  MuN2_2nd <- c()
  for(l in (G1+1):G){
    tmp <- MuN1[l] + truncdist::rtrunc(n = 1, 
                                       spec = 'norm',
                                       mean = 0,
                                       sd = 1.5,
                                       a = 0.1, 
                                       b = 3)^rbinom(1, size=1, prob=0.5)
    while(tmp <= 0) tmp <- MuN1[l] + truncdist::rtrunc(n = 1, 
                                                       spec = 'norm',
                                                       mean = 0,
                                                       sd = 1.5,
                                                       a = 0.1, 
                                                       b = 3)^rbinom(1, size=1, prob=0.5)
    MuN2_2nd <- c(MuN2_2nd, tmp)
  }
  MuN2 <- c(MuN2_1st, MuN2_2nd)
  MuT <- rnorm(G, 7, 1.5)
  Mu <- cbind(MuN1, MuN2, MuT)
  colnames(Mu) <- c('MuN1', 'MuN2', 'MuT')
  rownames(Mu) <- paste('Gene', seq = seq(1, G))
  
  ## Simulate SigmaN1, SigmaN2 and SigmaT for each gene
  SigmaN1 <- runif(n = G, min = 0.1, max = 0.8)
  SigmaN2 <- runif(n = G, min = 0.1, max = 0.8)
  SigmaT <- runif(n = G, min = 0.1, max = 0.8)
  Sigma <- cbind(SigmaN1, SigmaN2, SigmaT)
  colnames(Sigma) <- c('SigmaN1', 'SigmaN2', 'SigmaT')
  rownames(Sigma) <- paste('Gene', seq = seq(1, G))
  
  ## Initial values
  data.N1 <- matrix(0, G, M1)
  data.N2 <- matrix(0, G, M2)
  data.Y <- matrix(0, G, My)
  True.data.N1 <- matrix(0, G, My)
  True.data.N2 <- matrix(0, G, My)
  True.data.T <- matrix(0, G, My)
  ## Creat row and column name
  rownames(data.N1) <- rownames(data.N2) <-
    rownames(data.Y) <- rownames(True.data.N1) <-
    rownames(True.data.N2) <- rownames(True.data.T) <-
    paste('Gene', seq = seq(1, G))
  colnames(data.N1) <- paste('Sample', seq = seq(1, M1))
  colnames(data.N2) <- paste('Sample', seq = seq(1, M2))
  colnames(data.Y) <- colnames(True.data.N1) <-
    colnames(True.data.N2) <-
    colnames(True.data.T) <- paste('Sample', seq = seq(1, My))
  
  ## Simulate Tumor Proportion
  pi <- matrix(0, 3, My)
  pi[1,] <- runif(n = My, min = 0.01, max = 0.97)
  for(j in 1:My){
    pi[2, j] <- runif(n = 1, min = 0.01, max = 0.98 - pi[1,j])
    pi[3, j] <- 1 - sum(pi[,j])
  }
  rownames(pi) <- c('PiN1', 'PiN2', 'PiT')
  colnames(pi) <- paste('Sample', seq = seq(1, My))
  
  ## Simulate Data
  for(k in 1:G){
    
    data.N1[k,] <- 2^rnorm(M1, MuN1[k], SigmaN1[k]); # normal reference 1
    
    data.N2[k,] <- 2^rnorm(M2, MuN2[k], SigmaN2[k]); # normal reference 1
    
    True.data.T[k,] <- 2^rnorm(My, MuT[k], SigmaT[k]);  # True Tumor
    
    True.data.N1[k,] <- 2^rnorm(My, MuN1[k], SigmaN1[k]);  # True Normal 1
    
    True.data.N2[k,] <- 2^rnorm(My, MuN2[k], SigmaN2[k]);  # True Normal 1
    
    data.Y[k,] <- pi[1,]*True.data.N1[k,] + pi[2,]*True.data.N2[k,] +
      pi[3,]*True.data.T[k,] # Mixture Tumor
    
  }
  
  ## Transfer into bioconductor format
  data.Y <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(data.Y)),
                                 rowData = DataFrame(row.names = rownames(data.Y)),
                                 colData = DataFrame(row.names = colnames(data.Y)))
  data.N1 <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(data.N1)),
                                  rowData = DataFrame(row.names = rownames(data.N1)),
                                  colData = DataFrame(row.names = colnames(data.N1)))
  data.N2 <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(data.N2)),
                                  rowData = DataFrame(row.names = rownames(data.N2)),
                                  colData = DataFrame(row.names = colnames(data.N2)))
  True.data.T <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(True.data.T)),
                                      rowData = DataFrame(row.names = rownames(True.data.T)),
                                      colData = DataFrame(row.names = colnames(True.data.T)))
  True.data.N1 <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(True.data.N1)),
                                       rowData = DataFrame(row.names = rownames(True.data.N1)),
                                       colData = DataFrame(row.names = colnames(True.data.N1)))
  True.data.N2 <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(True.data.N2)),
                                       rowData = DataFrame(row.names = rownames(True.data.N2)),
                                       colData = DataFrame(row.names = colnames(True.data.N2)))
  
  if(output.more.info){
    test.data = list(pi=pi, Mu = Mu, Sigma = Sigma,
                     data.Y = data.Y, data.N1 = data.N1, data.N2 = data.N2,
                     True.data.N1 = True.data.N1, True.data.N2 = True.data.N2,
                     True.data.T = True.data.T)
  }else{
    test.data = list(pi=pi, Mu = Mu, Sigma = Sigma,
                     data.Y = data.Y, data.N1 = data.N1, data.N2 = data.N2)
  }
  
  
  return(test.data)
}
set.seed(123)
test.data.3comp = simulate_3comp(G1 = 675, G2 = 25, My = 20, M1 = 100, M2 = 100)
usethis::use_data(test.data.3comp, compress = 'xz', overwrite = TRUE)
