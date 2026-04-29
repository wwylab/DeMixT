#' @title Deconvolution of heterogeneous tumor samples with two components using
#' expression data from RNAseq or microarray platforms
#' 
#' @description DeMixNB is a variant of DeMixT such that the expressions are 
#' assumed to follow the negative binomial function 
#' 
#' @param data.Y A SummarizedExperiment object of expression data from mixed 
#' tumor samples. It is a \eqn{G} by \eqn{My} matrix where \eqn{G} is the number
#' of genes and \eqn{My} is the number of mixed samples. Samples with the same
#' tissue type should be placed together in columns.
#' 
#' @param data.N1 A SummarizedExperiment object of expression data 
#' from reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{M1} matrix 
#' where \eqn{G} is the number of genes and \eqn{M1} is the number of samples 
#' for component 1. 
#' 
#' @param niter The maximum number of iterations used in the algorithm of 
#' iterated conditional modes. A larger value better guarantees 
#' the convergence in estimation but increases the running time. The default is 
#' 30. 
#' 
#' @param nspikein The number of spikes in normal reference used for proportion
#' estimation. The default value is \eqn{ min(200, My)}, where 
#' \eqn{My} the number of mixed samples. If it is set to 0, proportion 
#' estimation is performed without any spike in normal reference.
#' 
#' @param tol The convergence criterion. The default is 10^(-8).
#' 
#' @param nthread The number of cores to use for parallelization. The default
#' is 1, i.e. no parallelization at the top level. Note that OpenMP is used 
#' under the hood when calculating the complete likelihood.
#' 
#' @param seeds A numeric vector of random seeds to use for each run. 
#' The DeMixNB algorithm will be run once for each seed. Default is c(123, 456, 789).
#' 
#' @param consistency_threshold The maximum allowed difference (range) 
#' between pi estimates across runs for a spot to be considered "consistent". 
#' Default is 0.1.
#' 
#' @param output.more.info A boolean flag controlling whether to return the 
#' estimations in each iteration. Default is False.
#'
#' @param data.scale The constant scale to be multiplied to the data. Default is 1.
#'
#' @return 
#' \item{pi_t_summary}{A data frame containing the estimated tumor proportions for
#' each run and a consistency flag at each sample indicating whether the runs agree 
#' with each other based on the maximum difference.}
#' 
#' \item{pi_n_summary}{A data frame containing the estimated non-tumor proportions 
#' for each run.}
#' 
#' \item{all_run_results}{A list containing all the detailed outputs of each run,
#' including the \eqn{mu} and \eqn{phi} estimates.}
#' 
#' @author Liyang Xie
#' @examples
#' 
#' # Example: estimate proportions for simulated two-component data 
#'   data(test.data.2comp)
#' # res.NB = DeMixNB(data.Y = test.data.2comp$data.Y, 
#' #                    data.N1 = test.data.2comp$data.N1,
#' #                    niter = 10, tol = 10^(-5))
#' 
#' @export

## Estimating function
DeMixNB <- function(data.Y, data.N1, niter = 30, nspikein = 200, tol = 1e-8,
                    nthread = 1, seeds = c(123, 456, 789), output.more.info = FALSE,
                    consistency_threshold = 0.1, data.scale = 1.0){
  ##----------------------------------------
  ## Step 1: Set up system and prepare input data
  ##----------------------------------------
  pbapply::pboptions(type="timer")
  ## Load C helpers for likelihood calculation
  #if(Sys.info()["sysname"] == "Windows")
  #{
  #  dyn.load('src/DeMixNB.dll')
  #} else
  #{
  #  dyn.load('src/DeMixNB.so')
  #}
  #
  ## Transferring data format
  data.Y <- SummarizedExperiment::assays(data.Y)[[1]] * data.scale
  data.N1 <- SummarizedExperiment::assays(data.N1)[[1]] * data.scale
  
  ## Create gene names by index if not already present
  if(is.null(rownames(data.Y)))
  {
    rownames(data.Y) <- paste('Gene', seq = seq(1, nrow(data.Y)))
    rownames(data.N1) <- paste('Gene', seq = seq(1, nrow(data.N1)))
  }
  
  ## Create column names by index if not already present
  if (is.null(colnames(data.Y))) {
    colnames(data.Y) <- paste('Sample', seq = ncol(data.N1) + seq(1, ncol(data.Y)))
    colnames(data.N1) <- paste('Sample', seq = seq(1, ncol(data.N1)))
  }
  gene.name <- rownames(data.Y)
  
  data.Y <- data.Y[gene.name,]
  data.N1 <- data.N1[gene.name,]
  
  if (!is.matrix(data.Y)) stop("data.Y must be a matrix.")
  if (!is.matrix(data.N1)) stop("data.N1 must be a matrix.")
  if (nrow(data.Y) != nrow(data.N1)) stop("data.Y and data.N1 must have the same number of rows (genes).")
  
  ##----------------------------------------
  ## Step 2: Run DeMixNB using different seeds and save the results
  ##----------------------------------------
  
  num_runs <- length(seeds)
  results_list <- vector("list", num_runs) # Pre-allocate list
  all_final_pi <- vector("list", num_runs)
  
  cat("Running DeMixNB", num_runs, "times with different seeds...\n")
  for (r in 1:num_runs) {
    current_seed <- seeds[r]
    cat("--- Starting Run", r, "/", num_runs, "with seed", current_seed, "---\n")
    set.seed(current_seed)
    
    # Run the core algorithm
    run_result <- tryCatch({
      DeMixNB_single_seed(data.Y = data.Y,
                          data.N1 = data.N1,
                          niter = niter,
                          nspikein = nspikein,
                          tol = tol,
                          nthread = nthread,
                          seed = current_seed,
                          output.more.info = output.more.info)
    }, error = function(e) {
      warning("Error occurred during run ", r, " with seed ", current_seed, ": ", e$message, call. = FALSE)
      NULL # Return NULL if an error occurs
    })
    
    results_list[[r]] <- run_result # Store result (or NULL)
    
    if (is.null(run_result)) {
      warning("Run ", r, " failed. Storing NA for pi estimates.")
      n_spots <- ncol(data.Y)
      all_final_pi[[r]] <- rep(NA_real_, n_spots)
    } else {
      all_final_pi[[r]] <- run_result$pi
    }
    
    cat("--- Finished Run", r, "---\n")
  }
  
  ##----------------------------------------
  ## Step 3: Check consistency of pi_t estimates and return results
  ##----------------------------------------
  cat("Combining results and checking consistency...\n")
  # Combine final pi estimates into a matrix (spots x runs)
  all_final_pi_t <- lapply(all_final_pi, function(mat) mat[2, ])
  all_final_pi_n <- lapply(all_final_pi, function(mat) mat[1, ])
  pi_t_matrix <- do.call(cbind, all_final_pi_t)
  pi_n_matrix <- do.call(cbind, all_final_pi_n)
  if (is.null(pi_t_matrix) || ncol(pi_t_matrix) == 0) {
    warning("No successful runs completed. Cannot create pi summary.")
    return(list(pi_summary = NULL, mu_phi_run1 = NULL, all_run_results = results_list))
  }
  # Ensure column names are set even if some runs failed and produced NAs
  run_indices <- 1:num_runs
  colnames(pi_t_matrix) <- paste0("pi_run", run_indices)
  colnames(pi_n_matrix) <- paste0("pi_run", run_indices)
  
  # Calculate max absolute difference (range) for each spot across runs, ignoring NAs
  pi_t_range <- apply(pi_t_matrix, 1, function(x) {
    finite_x <- x[is.finite(x)]
    if (length(finite_x) < 2) return(NA_real_) # Need at least 2 valid runs to calc range
    max(finite_x) - min(finite_x)
  })
  
  # Determine consistency (NA if range is NA or threshold is NA)
  consistency <- ifelse(is.na(pi_t_range) | is.na(consistency_threshold), NA_character_,
                        ifelse(pi_t_range > consistency_threshold, "inconsistent", "consistent"))
  
  # Create the final summary table
  pi_t_summary_table <- as.data.frame(pi_t_matrix)
  pi_t_summary_table$max_diff <- pi_t_range
  pi_t_summary_table$consistency <- consistency
  pi_n_summary_table <- as.data.frame(pi_n_matrix)
  
  ## Unload C helpers
  #if(Sys.info()["sysname"] == "Windows")
  #{
  #  dyn.unload('src/DeMixNB.dll')
  #} else
  #{
  #  dyn.unload('src/DeMixNB.so')
  #}
  #
  # Return the summary table and optionally the list of full results
  return(list(pi_t_summary = pi_t_summary_table,
              pi_n_summary = pi_n_summary_table,
              all_run_results = results_list)) 
}



#' @title Deconvolution of heterogeneous tumor samples with two components using
#' expression data from RNAseq or microarray platforms by single seed

#' @description The DeMixNB algorithm with a single seed.
#' @param data.Y A SummarizedExperiment object of expression data from mixed 
#' tumor samples. It is a \eqn{G} by \eqn{My} matrix where \eqn{G} is the number
#' of genes and \eqn{My} is the number of mixed samples. Samples with the same
#' tissue type should be placed together in columns.
#' 
#' @param data.N1 A SummarizedExperiment object of expression data 
#' from reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{M1} matrix 
#' where \eqn{G} is the number of genes and \eqn{M1} is the number of samples 
#' for component 1. 
#' 
#' @param niter The maximum number of iterations used in the algorithm of 
#' iterated conditional modes. A larger value better guarantees 
#' the convergence in estimation but increases the running time. The default is 
#' 30. 
#' 
#' @param nspikein The number of spikes in normal reference used for proportion
#' estimation. The default value is \eqn{ min(200, My)}, where 
#' \eqn{My} the number of mixed samples. If it is set to 0, proportion 
#' estimation is performed without any spike in normal reference.
#' 
#' @param tol The convergence criterion. The default is 10^(-8).
#' 
#' @param nthread The number of cores to use for parallelization. The default
#' is 1, i.e. no parallelization at the top level. Note that OpenMP is used 
#' under the hood when calculating the complete likelihood.
#' 
#' @param seed The random seed to use for each run. 
#' 
#' @param output.more.info A boolean flag controlling whether to return the 
#' estimations in each iteration. Default is False.
#' 
#' @return 
#' \item{pi}{A matrix of estimated proportion. First row corresponds to the 
#' proportion estimate for the known components and unknown component 
#' respectively, and each column 
#' corresponds to one sample.}
#' 
#' \item{pi.iter}{Estimated proportions in each iteration. It is a \eqn{niter*
#' My*2} array. This is enabled only when output.more.info = TRUE.}
#' 
#' \item{Mu}{A matrix of estimated \eqn{mu} of the negative binomial 
#' distribution for both components. Each row corresponds to one gene.} 
#' \item{Phi}{A matrix of estimated \eqn{phi} of the negative binomial 
#' distribution for both components. Each row corresponds to one gene.}
#' \item{gene.name}{The names of genes used in estimating the proportions. 
#' If no gene names are provided in the original data set, the genes will be
#' automatically indexed.}
DeMixNB_single_seed <- function(data.Y, data.N1, niter = 30, nspikein = 200, tol = 10e-8,
                                nthread = 1, output.more.info = FALSE, seed=100)
{
  # Set random seed
  set.seed(seed)
  # G: number of genes (rows in data.Y)
  # My: number of spots (columns in data.Y)
  G <- nrow(data.Y)
  My <- ncol(data.Y)
  gene.name <- rownames(data.Y)
  
  ##----------------------------------------
  ## Step 1: Estimate Normal Parameters (mu_n and phi_n) - Original Logic
  ##----------------------------------------
  cat("Step 1: Estimating Normal parameters (mu_n, phi_n)...\n") # Added status message
  mu_n <- phi_n <- rep(0, G) 
  for(g in 1:G){
    Fit_try <- try(fitdist(data.N1[g, ], distr = 'nbinom'), silent = TRUE)
    if("try-error" %in% class(Fit_try)){
      un <- mean(data.N1[g, ]) 
      vn <- var(data.N1[g, ])
      mu_n[g] <- un
      phi_n[g] <- ifelse(vn <= un | !is.finite(vn) | !is.finite(un), 1, un^2/(vn - un)) # Added checks for NA/NaN/Inf
      #if (phi_n[g] <= 0) phi_n[g] <- 1 # Ensure size is positive
    } else {
      Fit_NB <- fitdist(data.N1[g, ], distr = 'nbinom') 
      mu_n[g] <- Fit_NB$estimate[2] 
      phi_n[g] <- Fit_NB$estimate[1] 
    }
  }
  
  ##----------------------------------------
  ## Step 2 - Create spikein samples based on estimations of mu_n and phi_n
  ##----------------------------------------
  cat("Step 2: Generating spikein samples if needed\n")
  ## Generate Spikein samples
  nspikein <- min(nspikein, My)
  if (nspikein > 0)
  {
    Spikein.Normal = array(0, c(nrow(data.Y), nspikein))
    for (k in 1:nrow(data.Y))
    {
      Spikein.Normal[k,] <- rnbinom(nspikein, size = phi_n[k], mu = mu_n[k])
    }
    data.Y <- cbind(data.Y, Spikein.Normal)
  }
  Ny = ncol(data.Y)
  
  ##----------------------------------------
  ## Step 3: Initialize Tumor Parameters (pi, mu_t, phi_t) - Original Logic
  ##----------------------------------------
  cat("Step 3: Initializing Tumor parameters (pi_0, mu_t, phi_t)...\n") # Added status message
  # Original pi initialization
  pi_0 <- 0.5 + runif(Ny, min = -0.2, max = 0.2)
  
  mu_t <- phi_t <- rep(0, G)
  for(g in 1:G){
    ut <- mean(data.Y[g, ])
    vt <- var(data.Y[g, ])
    mu_t[g] <- ut
    phi_t[g] <- ifelse(vt <= ut | !is.finite(vt) | !is.finite(ut), 1, ut^2/(vt - ut))
  }
  
  ##----------------------------------------
  ## Step 4: Pre-allocate Storage for Iterative Results (Reduced)
  ##----------------------------------------
  cat("Step 4: Initializing result vectors \n")
  pi_iters    <- array(0, c(Ny, niter + 1)) # Stores pi estimates computed *in* iteration t (stored at t+1)
  mu_t_iters  <- array(0, c(G, niter + 1)) # Stores mu_t estimates computed *in* iteration t (stored at t+1)
  phi_t_iters <- array(0, c(G, niter + 1)) # Stores phi_t estimates computed *in* iteration t (stored at t+1)
  
  # Store initial estimates (iteration 0, these are used to start iteration 1)
  pi_iters[, 1]   <- pi_0
  mu_t_iters[, 1] <- mu_t
  phi_t_iters[, 1]<- phi_t
  
  ## Default: all iterations executed (will be truncated if convergence is reached)
  iters_used <- niter + 1
  
  ##----------------------------------------
  ## Step 5: Iterative Estimation Loop (Original Order/Logic)
  ##----------------------------------------
  cat("Step 5: Starting iterative estimation...\n") # Added status message
  for(t in 1:niter){
    cat('Iterations ', t, "/",  niter, '...\n')
    
    # Parameters from the START of this iteration (or initialization for t=1)
    pi_start_iter <- pi_iters[, t]
    mu_t_start_iter <- mu_t_iters[, t]
    phi_t_start_iter <- phi_t_iters[, t]
    
    ## (A) Estimate gene-specific parameters using pi from start of iteration.
    #cat("    Estimating gene parameters (mu_t, phi_t)...\n") # Added status message
    
    mu_phi_t_list <- pbapply::pblapply(X=1:G, FUN = Estimate_T_g, cl = nthread, 
                              data.Y = data.Y, pi_t = pi_start_iter, n = Ny,
                              mu_t = mu_t_start_iter, phi_t = phi_t_start_iter, # Pass starting estimates
                              mu_n = mu_n, phi_n = phi_n)
    # Check if the result is a list of vectors before binding
    is_list_of_vectors <- all(sapply(mu_phi_t_list, is.vector)) && all(sapply(mu_phi_t_list, length) == 2)
    if (!is_list_of_vectors || length(mu_phi_t_list) != G) {
      stop("Error in Estimate_T_g results: Expected a list of G vectors of length 2.")
    }
    mu_phi_t_matrix <- do.call(rbind, mu_phi_t_list)
    mu_t_new <- mu_phi_t_matrix[, 1]
    phi_t_new <- mu_phi_t_matrix[, 2]
    
    
    ## (B) Estimate tumor proportions pi using parameters from START of iteration
    #cat("    Estimating spot proportions (pi)...\n") # Added status message
    pi_t_list <- pbapply::pblapply(1:Ny, FUN = Estimate_pi_i, cl = nthread, 
                          data.Y = data.Y, pi_0 = pi_start_iter, G = G, 
                          mu_t = mu_t_start_iter, phi_t = phi_t_start_iter,
                          mu_n = mu_n, phi_n = phi_n)
    pi_t_new <- unlist(pi_t_list)
    pi_t_new <- pmax(0.05, pmin(pi_t_new, 0.95)) # ensure pi stays within [0.05, 0.95]
    
    
    # Store the newly computed parameters for this iteration t (in column t+1)
    pi_iters[, t + 1]   <- pi_t_new
    mu_t_iters[, t + 1] <- mu_t_new
    phi_t_iters[, t + 1]<- phi_t_new
    
    ## (D) Check convergence on pi: Compare pi at start of iteration (pi_start_iter) with new pi (pi_t_new)
    pi_change <- max(abs(pi_start_iter - pi_t_new))
    #cat("    Max pi change:", pi_change, "\n") # Added status message
    if(pi_change < tol){
      cat("Convergence reached at iteration", t, "\n") # Added status message
      iters_used <- t + 1 # Record the number of iterations used (column index)
      break
    }
  }
  
  if (iters_used == niter + 1 && niter > 0) {
    cat("Warning: Algorithm did not converge within ", niter, " iterations.\n")
  }
  
  ##----------------------------------------
  ## Step 6: Clean up and return results
  ##----------------------------------------
  cat("Step 6: Finalizing results...\n") # Added status message
  pi_iters    <- pi_iters[, 1:iters_used, drop = FALSE]
  mu_t_iters  <- mu_t_iters[, 1:iters_used, drop = FALSE]
  phi_t_iters <- phi_t_iters[, 1:iters_used, drop = FALSE]
  
  # Get the estimates from the *last computed* iteration (index iters_used)
  final_pi <- pi_iters[, iters_used]
  final_mu_t <- mu_t_iters[, iters_used]
  final_phi_t <- phi_t_iters[, iters_used]
  
  # Remove the pi estimates for the spikein samples
  pi_t <- final_pi[1:My]
  pi <- rbind(1 - pi_t, pi_t)
  colnames(pi) <- colnames(data.Y)[1:My]
  rownames(pi) <- c("pi_n", "pi_t")
  
  Mu <- cbind(mu_n, final_mu_t)
  rownames(Mu) <- gene.name
  
  phi <- cbind(phi_n, final_phi_t)
  rownames(phi) <- gene.name
  
  if (output.more.info)
  {
    return(list(pi = pi, pi_iters = pi_iters, 
                Mu = Mu, phi = phi,
                gene.name = gene.name))
  }
  else
  {
    return(list(pi = pi, 
                Mu = Mu, phi = phi,
                gene.name = gene.name))
  }
}

## Likelihood for pi_i
Likelihood_NB_pi_i <- function(pi_i, y_i, mu_t, phi_t, mu_n, phi_n, nG){
  
  LLk_pi <- 0
  res <- .C("LLK_NB_pi_i", as.double(pi_i), as.double(y_i), as.double(mu_t),
            as.double(phi_t), as.double(mu_n), as.double(phi_n), 
            as.integer(nG), 0)[[8]]
  return(res)
  
}

## Estimate pi_i
Estimate_pi_i <- function(i, data.Y, pi_0, G, mu_t, phi_t, mu_n, phi_n){
  
  y_i = data.Y[,i]
  
  Fit_opt <- optim(pi_0[i], fn = Likelihood_NB_pi_i,
                   method = 'L-BFGS-B',
                   lower = 0.05, upper = 0.95,
                   y_i = y_i, nG = G,
                   mu_t = mu_t, phi_t = phi_t,
                   mu_n = mu_n, phi_n = phi_n,
                   hessian = T
  )
  
  return(Fit_opt$par)
}

## Likelihood for gene-specific parameters
Likelihood_NB_T_g <- function(paras, y_g, pi, mu_ng, phi_ng, n = n){
  
  mu_tg = paras[1]; phi_tg = paras[2]
  LLk_T <- 1e-8
  res = .C("LLK_NB_T_g", as.double(paras),
           as.double(y_g), as.double(pi),
           as.double(mu_ng), as.double(phi_ng), as.integer(n), 0)[[7]]
  return(res)
}

## Estimate gene-specific para
Estimate_T_g <- function(g, data.Y, pi_t, n, mu_t, phi_t, mu_n, phi_n){
  
  y_g = data.Y[g,]
  Fit_opt_T <- optim(par = c(mu_t[g], phi_t[g]),
                     fn = Likelihood_NB_T_g,
                     method = 'L-BFGS-B',
                     lower = c(0.01, 0.01),
                     upper = c(1e8, 1e8),
                     y_g = y_g, n = n, pi = pi_t,
                     mu_ng = mu_n[g], phi_ng = phi_n[g],
                     hessian = T)
  mu_tg <- Fit_opt_T$par[1]
  phi_tg <- Fit_opt_T$par[2]
  return(c(mu_tg, phi_tg))
}

