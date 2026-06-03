# Load necessary libraries
library(MASS)
library(mvtnorm)
library(Matrix)
library(tictoc)
library(RSpectra)
library(RcppHungarian)
library(igraph)
library(flare)
library(polycor)
library(mclust)

library(foreach)
library(doParallel)

perm <- function(x, p, K) {
  x_perm <- rep(NA, length(x))
  for (i in 1:length(x)) {
    for (k in 1:K) {
      x_perm[x == k] <- p[k]
    }
  }
  return(x_perm)
}

# flatten ordinal data into a binary matrix 
flatten_ordinal <- function(R_mat, C_vec) {
  N <- nrow(R_mat)
  J <- ncol(R_mat)
  total_C <- sum(C_vec)
  R_flat <- matrix(0, nrow = N, ncol = total_C)
  
  col_start <- c(1, cumsum(C_vec)[-J] + 1)
  
  for (j in 1:J) {
    idx <- col_start[j] + R_mat[, j] - 1
    R_flat[cbind(1:N, idx)] <- 1
  }
  return(R_flat)
}

set.seed(123)

# ==============================================================================
# Set Model Parameters
# ==============================================================================
s <- 0.6  # edge probability: 0.2 or 0.6
N <- 1000 # Fixed sample size
J <- 100  # Total number of variables: 100 or 500
K <- 5    # Number of latent classes: 5 or 10
L <- 5    # Number of blocks: 5 or 10

block_sizes <- rep(J/L, L) 
V <- rep(seq_along(block_sizes), times = block_sizes)

pi_vec <- rep(1/K, K) # latent proportions
C_vec <- c(rep(2, J/2), rep(4, J/2)) # number of categories


N_rep <- 100

# ==============================================================================
# Main Simulation
# ==============================================================================
num_cores <- parallel::detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)

results_matrix <- foreach(n_rep = 1:N_rep, 
                          .combine = rbind, 
                          .packages = c("MASS", "mvtnorm", "Matrix", "tictoc", "RSpectra", "RcppHungarian", "igraph", "flare", "polycor", "huge", "mclust")) %dopar% {
                            
  tic()
  
  Sigma_list <- list()
  Omega_list <- list()
  chol_Sigma_list <- list() # Store sparse Cholesky factors for fast generation
  
  # Define Covariance/Precision Matrices
  for (k in 1:K) {
    prec_blocks <- list()
    blocks <- list()
    
    for (ell in 1:length(block_sizes)) {
      m <- block_sizes[ell]
      
      # Generate Random Sparse Graph
      B <- matrix(0, nrow = m, ncol = m)
      upper_idx <- which(upper.tri(B))
      
      B[upper_idx] <- rbinom(length(upper_idx), size = 1, prob = s) * 0.5
      B <- B + t(B)
      
      if (sum(abs(B)) == 0) {
        prec_block_mat <- diag(1, m)
      } else {
        ev <- eigen(B, symmetric = TRUE, only.values = TRUE)$values
        mu_max <- max(ev)
        mu_min <- min(ev)
        
        delta <- (mu_max - m * mu_min) / (m - 1)
        prec_block_mat <- B + diag(delta, m)
        
        D_inv <- diag(1 / sqrt(diag(prec_block_mat)))
        prec_block_mat <- D_inv %*% prec_block_mat %*% D_inv
      }
      
      prec_block_mat <- (prec_block_mat + t(prec_block_mat)) / 2
      
      # Final Transformation & Scaling
      block_mat_raw <- solve(prec_block_mat)
      block_mat_raw <- (block_mat_raw + t(block_mat_raw)) / 2
      
      block_mat <- cov2cor(block_mat_raw)
      block_mat <- (block_mat + t(block_mat)) / 2
      
      lambda_min <- min(eigen(block_mat, symmetric = TRUE, only.values = TRUE)$values)
      if (lambda_min < 1e-6) {
        block_mat <- block_mat + diag(abs(lambda_min) + 1e-4, m)
        block_mat <- cov2cor(block_mat)
        block_mat <- (block_mat + t(block_mat)) / 2
      }
      
      prec_block <- solve(block_mat)
      prec_block[B == 0 & row(prec_block) != col(prec_block)] <- 0
      
      prec_blocks[[ell]] <- prec_block
      blocks[[ell]] <- block_mat
    }
    
    Omega_list[[k]] <- as.matrix(bdiag(prec_blocks))
    Sigma_list[[k]] <- as.matrix(bdiag(blocks))
    
    chol_blocks <- lapply(blocks, chol)
    chol_Sigma_list[[k]] <- bdiag(chol_blocks)
  }
  
  # Define Thresholds and Generate Data
  mu_list <- list()          # thresholds for binary variables
  thresh_poly_list <- list() # thresholds for polytomous variables
  
  for (k in 1:K) {
    mu_list[[k]] <- runif(J/2, -1, 1)
    thresh_poly_list[[k]] <- t(replicate(J/2, sort(runif(3, -2, 2))))
  }
  
  X_mat <- matrix(NA, nrow = N, ncol = J)
  R_mat <- matrix(NA, nrow = N, ncol = J)
  
  # Latent Cluster Z
  Z_vec <- sample(1:K, size = N, replace = TRUE, prob = pi_vec)
  
  for (i in 1:N) {
    k <- Z_vec[i]
    
    # Latent Gaussian X
    z_i <- rnorm(J)
    x_i <- as.numeric(z_i %*% chol_Sigma_list[[k]])
    X_mat[i, ] <- x_i
    
    # ordinal responses R
    r_i_bin <- ifelse(x_i[1:(J/2)] > mu_list[[k]], 2, 1)
    thresh_poly <- thresh_poly_list[[k]]
    r_i_poly <- sapply(1:(J/2), function(j) {
      findInterval(x_i[J/2 + j], vec = c(-Inf, thresh_poly[j, ], Inf))
    })
    
    R_mat[i, ] <- c(r_i_bin, r_i_poly)
  }
  
  # ==============================================================================
  # Step 1: Spectral Clustering with Flattened Ordinal Data
  # ==============================================================================
  t0_step1 <- Sys.time()
  
  R_flat <- flatten_ordinal(R_mat, C_vec)
  
  nstart <- 50
  U <- svds(R_flat, k=K)$u
  kmeans_res <- kmeans(U, centers=K, iter.max=100, nstart=nstart)
  
  Z_hat <- kmeans_res$cluster
  
  # Map estimated labels to true labels
  mat <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      mat[i, j] <- sum((Z_hat == i) & (Z_vec == j))
    }
  }
  perm_mat <- HungarianSolver(-mat)$pairs[, 2]
  Z_hat <- perm(Z_hat, perm_mat, K)
  
  rep_acc_Z <- sum(Z_hat != Z_vec)
  
  time_step1 <- as.numeric(difftime(Sys.time(), t0_step1, units = "secs"))
  
  # ==============================================================================
  # Step 2: Covariance Matrix & Block Structure Estimation
  # ==============================================================================
  t0_step2 <- Sys.time()
  
  Sigma_hat_list <- list()
  lambda_cov <- sqrt(log(J*K)/N)
  
  for (k in 1:K) {
    indices_k <- which(Z_hat == k)
    R_k <- R_mat[indices_k, , drop = FALSE]
    N_k <- length(indices_k)
    
    if (N_k < 2) {
      Sigma_hat_list[[k]] <- matrix(0, J, J)
      next
    }
    
    pc_cor <- tryCatch({
      hetcor(R_k, ML = FALSE, std.err = FALSE)$correlations
    }, error = function(e) cor(R_k))
    
    Sigma_k <- pc_cor
    Sigma_k[is.na(Sigma_k)] <- 0
    
    # Thresholding (if we do indicator instead)
    # Sigma_k[abs(Sigma_k) < lambda_cov] <- 0 ### this helps for large L
    # diag(Sigma_k) <- 1
    Sigma_hat_list[[k]] <- Sigma_k
  }
  
  # Estimate Block Structures via Leiden
  sum_I <- matrix(0, nrow = J, ncol = J)
  for (k in 1:K) {
    sum_I <- sum_I + abs(Sigma_hat_list[[k]]) 
  }
  diag(sum_I) <- 0
  avg_support <- sum_I / K
  
  g <- graph_from_adjacency_matrix(avg_support, mode = "undirected", weighted = TRUE)
  leiden_res <- cluster_leiden(g, objective_function = "modularity")
  
  L_hat <- length(leiden_res)
  V_hat <- membership(leiden_res)
  
  rep_acc_L <- as.numeric(L_hat != L)
  rep_acc_V <- adjustedRandIndex(V, V_hat)
  time_step2 <- as.numeric(difftime(Sys.time(), t0_step2, units = "secs"))

  
  # ==============================================================================
  # Step 3: Precision Matrix Estimation and tuning selection via CV
  # ==============================================================================
  t0_step3 <- Sys.time()

  lambda_vec <- c(2, 1, 0.5, 0.1) * sqrt(log(J*K)/N) 
  n_lambda <- length(lambda_vec)
  L_labels <- sort(unique(V_hat))
  
  D_folds <- 5
  folds <- sample(rep(1:D_folds, length.out = N))
  
  cv_loss <- rep(0, n_lambda)
  
  for (d in 1:D_folds) {
    train_idx <- which(folds != d)
    test_idx <- which(folds == d)
    
    for (k in 1:K) {
      train_k <- intersect(train_idx, which(Z_hat == k))
      test_k <- intersect(test_idx, which(Z_hat == k))
      N_test_k <- length(test_k)
      
      # estimate train covariance matrix
      R_train_k <- R_mat[train_k, , drop = FALSE]
      Sigma_train_k <- hetcor(R_train_k, ML = FALSE, std.err = FALSE)$correlations
      Sigma_train_k[is.na(Sigma_train_k)] <- 0
      Sigma_train_k[abs(Sigma_train_k) < lambda_cov] <- 0 # Apply threshold
      diag(Sigma_train_k) <- 1
      
      # estimate test covariance matrix
      R_test_k <- R_mat[test_k, , drop = FALSE]
      Sigma_test_k <- hetcor(R_test_k, ML = FALSE, std.err = FALSE)$correlations
      Sigma_test_k[is.na(Sigma_test_k)] <- 0
      diag(Sigma_test_k) <- 1
      
      # precision matrix estimation
      for (l_idx in seq_along(L_labels)) {
        l <- L_labels[l_idx]
        indices <- which(V_hat == l)
        
        if (length(indices) == 0) next
        
        Sigma_sub_train <- Sigma_train_k[indices, indices, drop = FALSE]
        Sigma_sub_test <- Sigma_test_k[indices, indices, drop = FALSE]
        
        # Ensure PD for the training sub-matrix
        Sigma_sub_train <- (Sigma_sub_train + t(Sigma_sub_train)) / 2
        min_eig <- min(eigen(Sigma_sub_train, symmetric = TRUE, only.values = TRUE)$values)
        if (min_eig < 1e-4) {
          Sigma_sub_train <- Sigma_sub_train + diag(abs(min_eig) + 1e-4, nrow(Sigma_sub_train))
        }
        
        # Fit GLASSO on the training block
        out <- huge(Sigma_sub_train, method = "glasso", lambda = lambda_vec, verbose = FALSE)
        
        # Evaluate the loss
        for (lambda_ind in 1:n_lambda) {
          Omega_kl_hat <- as.matrix(out$icov[[lambda_ind]])
          tr_val <- sum(diag(Sigma_sub_test %*% Omega_kl_hat))

          ev <- eigen(Omega_kl_hat, symmetric = TRUE, only.values = TRUE)$values
          logdet <- if (any(ev <= 1e-10)) sum(log(pmax(ev, 1e-10))) else sum(log(ev))
          
          cv_loss[lambda_ind] <- cv_loss[lambda_ind] + N_test_k * (tr_val - logdet)
        }
      }
    }
  }
  
  # Select optimal lambda
  opt_lambda_ind <- which.min(cv_loss)
  opt_lambda <- lambda_vec[opt_lambda_ind]
  
  time_step3 <- as.numeric(difftime(Sys.time(), t0_step3, units = "secs"))
  
  curr_timer <- toc(log = TRUE, quiet = TRUE)
  time_elapsed <- curr_timer$toc - curr_timer$tic
  
  # Return accuracy metrics and selected lambda
  c(acc_Z = rep_acc_Z, 
    acc_L = rep_acc_L,
    L_hat = L_hat,
    acc_V = rep_acc_V,
    time_step1 = time_step1,
    time_step2 = time_step2,
    time_step3 = time_step3,
    time_total = time_elapsed,
    lambda_ind = opt_lambda_ind
  )
}

stopCluster(cl)

# ==============================================================================
# Save Results
# ==============================================================================
results_df <- as.data.frame(results_matrix)

cat("\n=== Final Averaged Results ===\n")
cat("Misclustering Z (Hamming): ", mean(results_df$acc_Z), "\n")
cat("Error in L_hat:            ", mean(results_df$acc_L), "\n")
cat("Error in L_hat:            ", mean(results_df$L_hat) - L, "\n")
cat("1-ARI for Blocks (V_hat):    ", 1-mean(results_df$acc_V), "\n")

ind <- unique(results_df$lambda_ind)
ind_select = ind[which.max(tabulate(match(results_df$lambda_ind, ind)))]

final_summary <- data.frame(
  N = N,
  J = J,
  K = K,
  L = L,
  s = s,
  ind_select = ind_select
)
