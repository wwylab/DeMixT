#' @title Function to simulate two-component test data
#' @description Function to simulate two-component test data for DeMixT.
#'
#' @param G Number of genes for simulation.
#' @param My Number of mixture tumor samples for simulation.
#' @param M1 Number of normal reference for simulation.
#' @param output.more.info The logical flag indicating wheter to show True.data.T
#' and True.data.N1 in the output. The default is FALSE.
#'
#' @return
#' \item{pi}{A matrix of estimated proportion. First row and second row 
#' corresponds to the proportion estimate for the known components and unkown 
#' component respectively for two or three component settings. Each column 
#' corresponds to one sample.}
#' \item{Mu}{Simulated \eqn{Mu} of log2-normal distribution for both known
#' (\eqn{MuN1}) and unknown component (\eqn{MuT}).}
#' \item{Sigma}{Simulated \eqn{Sigma} of log2-normal distribution for both 
#' known (\eqn{SigmaN1}) and unknown component (\eqn{SigmaT}).}
#' \item{data.Y}{A SummarizedExperiment object of expression data from mixed 
#' tumor samples. It is a \eqn{G} by \eqn{My} matrix where \eqn{G} is the number
#' of genes and \eqn{My} is the number of mixed samples. Samples with the same
#' tissue type should be placed together in columns.}
#' \item{data.N1}{A SummarizedExperiment object of expression data 
#' from reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{M1} matrix 
#' where \eqn{G} is the number of genes and \eqn{M1} is the number of samples 
#' for component 1.} 
#' \item{True.data.T}{A SummarizedExperiment object of simulated tumor expression 
#' data. It is a \eqn{G} by \eqn{My} matrix, where \eqn{G} is the number of 
#' genes and \eqn{My} is the number of mixed samples. This is enabled only when 
#' output.more.info = TRUE.}
#' \item{True.data.N1}{A SummarizedExperiment object of simulated true 
#' expression data for reference component 1 (e.g., normal). It is a \eqn{G} 
#' by \eqn{M1} matrix where \eqn{G} is the number of genes and \eqn{M1} is the 
#' number of samples for component 1. This is enabled only when 
#' output.more.info = TRUE.}
#'
#' @name simulate_2comp
#' 
#' 
#' @export
#'
#' @examples
#' test.data = simulate_2comp(G = 500, My = 100, M1 = 100)
#' test.data$pi
#' test.data$Mu
#' test.data$Sigma
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

