# Load local R library if this repository has one.
script_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                        error = function(e) NA_character_)
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(script_file) else getwd()
local_libs <- c(file.path(script_dir, "LCM_simulation", "Rlib"),
                file.path(script_dir, "Rlib"))
local_libs <- local_libs[dir.exists(local_libs)]
if (length(local_libs) > 0) {
  .libPaths(unique(c(local_libs, .libPaths())))
}

required_packages <- c(
  "Matrix",
  "RSpectra",
  "RcppHungarian",
  "igraph",
  "polycor",
  "glasso",
  "EstimateGroupNetwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall EstimateGroupNetwork with:\n",
    "install.packages('EstimateGroupNetwork', repos = c('https://cran.r-universe.dev', ",
    "'https://cloud.r-project.org'))",
    call. = FALSE
  )
}

as_bool <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

env_int <- function(name, default) {
  as.integer(Sys.getenv(name, unset = as.character(default)))
}

env_num <- function(name, default) {
  as.numeric(Sys.getenv(name, unset = as.character(default)))
}

# Helper function to permute cluster labels to match the true labels.
perm <- function(x, p, K) {
  x_perm <- rep(NA_integer_, length(x))
  for (k in seq_len(K)) {
    x_perm[x == k] <- p[k]
  }
  x_perm
}

# Helper function to flatten ordinal data into a binary matrix
# for spectral clustering (as per Eq 7 in the manuscript).
flatten_ordinal <- function(R_mat, C_vec) {
  N <- nrow(R_mat)
  J <- ncol(R_mat)
  total_C <- sum(C_vec)
  R_flat <- matrix(0, nrow = N, ncol = total_C)

  col_start <- c(1, cumsum(C_vec)[-J] + 1)

  for (j in seq_len(J)) {
    idx <- col_start[j] + R_mat[, j] - 1
    R_flat[cbind(seq_len(N), idx)] <- 1
  }
  R_flat
}

choose2 <- function(x) {
  x * (x - 1) / 2
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n < 2) return(1)

  sum_comb <- sum(choose2(tab))
  row_comb <- sum(choose2(rowSums(tab)))
  col_comb <- sum(choose2(colSums(tab)))
  total_comb <- choose2(n)
  expected <- row_comb * col_comb / total_comb
  max_index <- (row_comb + col_comb) / 2
  denom <- max_index - expected

  if (abs(denom) < .Machine$double.eps) {
    return(as.numeric(identical(as.integer(factor(x)), as.integer(factor(y)))))
  }

  (sum_comb - expected) / denom
}

matched_block_accuracy <- function(V_true, V_hat) {
  true_labels <- sort(unique(V_true))
  hat_labels <- sort(unique(V_hat))
  M_max <- max(length(true_labels), length(hat_labels))
  mat_v <- matrix(0, nrow = M_max, ncol = M_max)

  for (i in seq_along(hat_labels)) {
    for (j in seq_along(true_labels)) {
      mat_v[i, j] <- sum((V_hat == hat_labels[i]) & (V_true == true_labels[j]))
    }
  }

  match_res <- RcppHungarian::HungarianSolver(-mat_v)$pairs
  sum(mat_v[match_res]) / length(V_true)
}

block_metrics <- function(V_true, V_hat, L_true) {
  L_hat <- length(unique(V_hat))
  c(
    acc_L = as.numeric(L_hat != L_true),
    L_hat = L_hat,
    acc_V = adjusted_rand_index(V_true, V_hat),
    acc_V_hamming = matched_block_accuracy(V_true, V_hat)
  )
}

sanitize_covariance <- function(Sigma, min_eig_target = 1e-4, force_correlation = TRUE) {
  Sigma <- as.matrix(Sigma)
  p <- ncol(Sigma)
  Sigma[!is.finite(Sigma)] <- 0
  Sigma <- (Sigma + t(Sigma)) / 2

  if (force_correlation) {
    diag(Sigma) <- 1
  } else {
    diag(Sigma) <- pmax(diag(Sigma), min_eig_target)
  }

  eig_vals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig_vals)
  if (!is.finite(min_eig) || min_eig < min_eig_target) {
    Sigma <- Sigma + diag(abs(min_eig) + min_eig_target, p)
    if (force_correlation) {
      Sigma <- cov2cor(Sigma)
      Sigma[!is.finite(Sigma)] <- 0
      diag(Sigma) <- 1
    }
  }

  (Sigma + t(Sigma)) / 2
}

estimate_glasso_precision <- function(Sigma, rho, precision_tol = 1e-6) {
  Sigma <- sanitize_covariance(Sigma)
  p <- ncol(Sigma)

  if (p == 1) {
    return(matrix(1 / max(Sigma[1, 1], 1e-10), nrow = 1, ncol = 1))
  }

  fit <- tryCatch(
    glasso::glasso(
      s = Sigma,
      rho = rho,
      penalize.diagonal = FALSE,
      trace = FALSE
    ),
    error = function(e) {
      Sigma_retry <- sanitize_covariance(Sigma, min_eig_target = 1e-2)
      glasso::glasso(
        s = Sigma_retry,
        rho = rho,
        penalize.diagonal = FALSE,
        trace = FALSE
      )
    }
  )

  Theta <- as.matrix(fit$wi)
  Theta <- (Theta + t(Theta)) / 2
  Theta[abs(Theta) < precision_tol] <- 0
  Theta
}

estimate_leiden_blocks <- function(Sigma_hat_list) {
  J <- nrow(Sigma_hat_list[[1]])
  K <- length(Sigma_hat_list)
  sum_I <- matrix(0, nrow = J, ncol = J)

  for (k in seq_len(K)) {
    sum_I <- sum_I + abs(Sigma_hat_list[[k]])
  }

  diag(sum_I) <- 0
  avg_support <- sum_I / K

  if (sum(avg_support) <= 0) {
    return(seq_len(J))
  }

  g <- igraph::graph_from_adjacency_matrix(
    avg_support,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  leiden_res <- tryCatch(
    igraph::cluster_leiden(g, objective_function = "modularity", resolution = 1),
    error = function(e) igraph::components(g)
  )

  as.integer(igraph::membership(leiden_res))
}

estimate_blockwise_glasso <- function(Sigma_hat_list, V_hat, rho) {
  K <- length(Sigma_hat_list)
  J <- nrow(Sigma_hat_list[[1]])
  L_labels <- sort(unique(V_hat))
  Omega_hat <- vector("list", length = K)

  for (k in seq_len(K)) {
    Omega_hat[[k]] <- matrix(0, nrow = J, ncol = J)
    Sigma_k <- Sigma_hat_list[[k]]

    for (l in L_labels) {
      indices <- which(V_hat == l)
      if (length(indices) == 0) next

      Sigma_sub <- Sigma_k[indices, indices, drop = FALSE]
      Omega_hat[[k]][indices, indices] <- estimate_glasso_precision(Sigma_sub, rho)
    }
  }

  Omega_hat
}

extract_blocks_from_precision <- function(Omega_hat_list, support_tol = 1e-6) {
  J <- nrow(Omega_hat_list[[1]])
  support <- matrix(FALSE, nrow = J, ncol = J)

  for (Omega_hat in Omega_hat_list) {
    A <- abs(as.matrix(Omega_hat)) > support_tol
    diag(A) <- FALSE
    support <- support | A
  }

  g <- igraph::graph_from_adjacency_matrix(
    support * 1,
    mode = "undirected",
    diag = FALSE
  )
  as.integer(igraph::components(g)$membership)
}

extract_blocks_from_jgl_screening <- function(Sigma_hat_list, nvec,
                                              lambda1, lambda2,
                                              penalty, weights) {
  S <- lapply(Sigma_hat_list, sanitize_covariance)
  K <- length(S)
  p <- ncol(S[[1]])

  weights_vec <- switch(
    weights,
    equal = rep(1, K),
    sample.size = nvec / sum(nvec)
  )

  if (penalty == "fused") {
    if (K == 2) {
      critboth <- matrix(FALSE, p, p)
      S_sum <- matrix(0, p, p)

      for (k in seq_len(K)) {
        critboth <- critboth | (abs(S[[k]]) * weights_vec[k] > lambda1 + lambda2)
        S_sum <- S_sum + weights_vec[k] * S[[k]]
      }

      critboth <- critboth | (abs(S_sum) > 2 * lambda1)
    } else {
      critboth <- matrix(FALSE, p, p)

      for (k in seq_len(K)) {
        critboth <- critboth | (abs(S[[k]]) * weights_vec[k] > lambda1)
      }
    }
  } else if (penalty == "group") {
    tempsum <- matrix(0, p, p)

    for (k in seq_len(K)) {
      tempsum <- tempsum + (pmax(weights_vec[k] * abs(S[[k]]) - lambda1, 0))^2
    }

    critboth <- tempsum > lambda2^2
  } else {
    stop("Unknown JGL penalty: ", penalty, call. = FALSE)
  }

  diag(critboth) <- FALSE
  g <- igraph::graph_from_adjacency_matrix(
    critboth * 1,
    mode = "undirected",
    diag = FALSE
  )
  as.integer(igraph::components(g)$membership)
}

precision_metrics <- function(Omega_hat_list, Omega_true_list,
                              true_zero_tol = 0.01, est_zero_tol = 1e-6) {
  K <- length(Omega_true_list)
  J <- nrow(Omega_true_list[[1]])
  acc <- rep(NA_real_, K)
  acc_F <- rep(NA_real_, K)

  for (k in seq_len(K)) {
    Omega_hat <- as.matrix(Omega_hat_list[[k]])
    Omega_true <- as.matrix(Omega_true_list[[k]])
    Omega_true[abs(Omega_true) < true_zero_tol] <- 0
    Omega_hat[abs(Omega_hat) < est_zero_tol] <- 0

    acc[k] <- mean((Omega_true == 0) == (Omega_hat == 0))
    acc_F[k] <- norm(Omega_true - Omega_hat, "F")^2 / J^2
  }

  c(
    acc_Omega = 1 - mean(acc),
    acc_Omega_F = sqrt(mean(acc_F))
  )
}

fit_estimate_group_network <- function(Sigma_hat_list, nvec,
                                       fixed_lambda1, nlambda2,
                                       lambda2_min_ratio, logseql2,
                                       criterion, gamma, dec,
                                       penalty, weights,
                                       maxiter, truncate,
                                       support_tol) {
  Rlist <- lapply(Sigma_hat_list, sanitize_covariance)
  names(Rlist) <- paste0("class", seq_along(Rlist))
  nvec <- as.integer(nvec)
  K <- length(Rlist)
  weights_vec <- switch(
    weights,
    equal = rep(1, K),
    sample.size = nvec / sum(nvec)
  )

  tol <- 10^(-dec)
  jgl_ns <- asNamespace("EstimateGroupNetwork")
  l2sequence <- get("l2sequence", jgl_ns)
  myJGL <- get("myJGL", jgl_ns)
  information_criterion <- get("InformationCriterion", jgl_ns)

  myJGLarglist <- list(
    penalty = penalty,
    weights = weights_vec,
    penalize.diagonal = FALSE,
    maxiter = maxiter,
    tol = tol,
    rho = 1,
    truncate = truncate
  )

  lambda2_grid <- l2sequence(
    S = Rlist,
    n = nvec,
    l1cand = fixed_lambda1,
    nlambda2 = nlambda2,
    lambda2.min.ratio = lambda2_min_ratio,
    logseql2 = logseql2,
    l2max = 1,
    ncores = 1,
    myJGLarglist = myJGLarglist
  )
  lambda2_grid <- unique(as.numeric(lambda2_grid))
  lambda2_grid <- lambda2_grid[is.finite(lambda2_grid) & lambda2_grid >= 0]
  if (length(lambda2_grid) == 0) {
    stop("The JGL lambda2 grid is empty.", call. = FALSE)
  }

  fits <- vector("list", length(lambda2_grid))
  ic_values <- rep(Inf, length(lambda2_grid))
  for (i in seq_along(lambda2_grid)) {
    fits[[i]] <- tryCatch(
      myJGL(
        S = Rlist,
        n = nvec,
        lambda1 = fixed_lambda1,
        lambda2 = lambda2_grid[i],
        penalty = penalty,
        weights = weights_vec,
        penalize.diagonal = FALSE,
        maxiter = maxiter,
        tol = tol,
        rho = 1,
        truncate = truncate
      ),
      error = function(e) NULL
    )

    if (!is.null(fits[[i]])) {
      ic_values[i] <- information_criterion(
        theta = fits[[i]]$concentrationMatrix,
        S = Rlist,
        n = nvec,
        criterion = criterion,
        count.unique = FALSE,
        gamma = gamma,
        dec = dec
      )
    }
  }

  if (all(!is.finite(ic_values))) {
    stop("All JGL fits failed while selecting lambda2.", call. = FALSE)
  }

  best_idx <- which.min(ic_values)
  lambda1 <- fixed_lambda1
  lambda2 <- lambda2_grid[best_idx]
  fit_raw <- fits[[best_idx]]

  fit <- list(
    network = fit_raw$network,
    concentrationMatrix = fit_raw$concentrationMatrix,
    correlationMatrix = Rlist,
    InformationCriteria = setNames(ic_values, paste0("lambda2=", signif(lambda2_grid, 4))),
    TuningParameters = c(lambda1 = lambda1, lambda2 = lambda2),
    Lambda2Grid = lambda2_grid,
    n = nvec,
    miscellaneous = c(
      method = "InformationCriterion",
      `information criterion used` = criterion,
      `fixed lambda1` = fixed_lambda1,
      nlambda2 = nlambda2,
      lambda2.min.ratio = lambda2_min_ratio,
      logseql2 = logseql2,
      `type of penalty` = ifelse(penalty == "fused", "Fused Graphical Lasso",
                                 "Group Graphical Lasso"),
      weights = weights,
      maxiter = maxiter,
      truncate = truncate
    )
  )

  names(fit$network) <- names(fit$concentrationMatrix) <- names(Rlist)

  Omega_hat <- lapply(fit$concentrationMatrix, function(x) {
    x <- as.matrix(x)
    x <- (x + t(x)) / 2
    x[abs(x) < support_tol] <- 0
    x
  })

  V_hat <- extract_blocks_from_jgl_screening(
    Sigma_hat_list = Rlist,
    nvec = nvec,
    lambda1 = lambda1,
    lambda2 = lambda2,
    penalty = penalty,
    weights = weights
  )

  list(
    Omega_hat = Omega_hat,
    V_hat = V_hat,
    lambda1 = lambda1,
    lambda2 = lambda2,
    fit = fit
  )
}

generate_sparse_model <- function(K, L, J, s, block_sizes) {
  Sigma_list <- vector("list", length = K)
  Omega_list <- vector("list", length = K)
  chol_Sigma_list <- vector("list", length = K)

  for (k in seq_len(K)) {
    prec_blocks <- vector("list", length = L)
    blocks <- vector("list", length = L)

    for (ell in seq_along(block_sizes)) {
      m <- block_sizes[ell]
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

      block_mat_raw <- solve(prec_block_mat)
      block_mat_raw <- (block_mat_raw + t(block_mat_raw)) / 2

      block_mat <- cov2cor(block_mat_raw)
      block_mat <- sanitize_covariance(block_mat)

      prec_block <- solve(block_mat)
      prec_block[B == 0 & row(prec_block) != col(prec_block)] <- 0

      prec_blocks[[ell]] <- prec_block
      blocks[[ell]] <- block_mat
    }

    Omega_list[[k]] <- as.matrix(Matrix::bdiag(prec_blocks))
    Sigma_list[[k]] <- as.matrix(Matrix::bdiag(blocks))
    chol_Sigma_list[[k]] <- as.matrix(Matrix::bdiag(lapply(blocks, chol)))
  }

  list(
    Sigma_list = Sigma_list,
    Omega_list = Omega_list,
    chol_Sigma_list = chol_Sigma_list
  )
}

run_replication <- function(n_rep) {
  set.seed(base_seed + n_rep)
  t_total <- proc.time()[["elapsed"]]
  if (show_progress) {
    cat(sprintf("Starting replication %d/%d\n", n_rep, N_rep))
    flush.console()
  }

  generated <- generate_sparse_model(K, L, J, s, block_sizes)
  Omega_list <- generated$Omega_list
  chol_Sigma_list <- generated$chol_Sigma_list

  # Generate ordinal data from the latent Gaussian model.
  mu_list <- vector("list", length = K)
  thresh_poly_list <- vector("list", length = K)

  for (k in seq_len(K)) {
    mu_list[[k]] <- runif(J / 2, -C_mu, C_mu)
    thresh_poly_list[[k]] <- t(replicate(J / 2, sort(runif(3, -2, 2))))
  }

  X_mat <- matrix(NA_real_, nrow = N, ncol = J)
  R_mat <- matrix(NA_integer_, nrow = N, ncol = J)
  Z_vec <- sample(seq_len(K), size = N, replace = TRUE, prob = pi_vec)

  for (i in seq_len(N)) {
    k <- Z_vec[i]
    z_i <- rnorm(J)
    x_i <- as.numeric(z_i %*% chol_Sigma_list[[k]])
    X_mat[i, ] <- x_i

    r_i_bin <- ifelse(x_i[seq_len(J / 2)] > mu_list[[k]], 2, 1)
    thresh_poly <- thresh_poly_list[[k]]
    r_i_poly <- sapply(seq_len(J / 2), function(j) {
      findInterval(x_i[J / 2 + j], vec = c(-Inf, thresh_poly[j, ], Inf))
    })

    R_mat[i, ] <- c(r_i_bin, r_i_poly)
  }

  # Step 1: Spectral clustering.
  t0_step1 <- proc.time()[["elapsed"]]

  R_flat <- flatten_ordinal(R_mat, C_vec)
  U <- RSpectra::svds(R_flat, k = K)$u
  kmeans_res <- kmeans(U, centers = K, iter.max = 100, nstart = nstart)
  Z_hat <- kmeans_res$cluster

  mat <- matrix(0, K, K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      mat[i, j] <- sum((Z_hat == i) & (Z_vec == j))
    }
  }
  perm_mat <- RcppHungarian::HungarianSolver(-mat)$pairs[, 2]
  Z_hat <- perm(Z_hat, perm_mat, K)

  rep_acc_Z <- sum(Z_hat != Z_vec)
  time_step1 <- proc.time()[["elapsed"]] - t0_step1

  # Step 2: Estimate the K class covariance/correlation matrices.
  t0_step2 <- proc.time()[["elapsed"]]

  Sigma_hat_list <- vector("list", length = K)
  nvec <- integer(K)
  lambda_cov <- sqrt(log(J * K) / N)

  for (k in seq_len(K)) {
    indices_k <- which(Z_hat == k)
    R_k <- R_mat[indices_k, , drop = FALSE]
    N_k <- length(indices_k)
    nvec[k] <- max(N_k, 2L)

    if (N_k < 2) {
      Sigma_hat_list[[k]] <- diag(1, J)
      next
    }

    pc_cor <- tryCatch(
      polycor::hetcor(R_k, ML = FALSE, std.err = FALSE)$correlations,
      error = function(e) cor(R_k)
    )

    Sigma_k <- as.matrix(pc_cor)
    Sigma_k[!is.finite(Sigma_k)] <- 0
    Sigma_k <- (Sigma_k + t(Sigma_k)) / 2
    diag(Sigma_k) <- 1

    Sigma_k[abs(Sigma_k) < lambda_cov] <- 0
    diag(Sigma_k) <- 1
    Sigma_hat_list[[k]] <- sanitize_covariance(Sigma_k)
  }

  time_step2 <- proc.time()[["elapsed"]] - t0_step2

  # Proposed LCM method: Leiden block recovery, then Step 3 glasso.
  t0_ours_block <- proc.time()[["elapsed"]]
  V_hat_ours <- estimate_leiden_blocks(Sigma_hat_list)
  ours_block <- block_metrics(V, V_hat_ours, L)
  time_ours_block <- proc.time()[["elapsed"]] - t0_ours_block

  t0_ours_precision <- proc.time()[["elapsed"]]
  Omega_hat_ours <- estimate_blockwise_glasso(Sigma_hat_list, V_hat_ours, opt_lambda)
  ours_precision <- precision_metrics(Omega_hat_ours, Omega_list)
  time_ours_precision <- proc.time()[["elapsed"]] - t0_ours_precision

  # Comparison method: EstimateGroupNetwork/JGL directly on the K covariance
  # matrices from Step 2. Its precision support gives the block components.
  t0_jgl <- proc.time()[["elapsed"]]
  jgl_status <- 1
  jgl_block <- c(acc_L = NA_real_, L_hat = NA_real_, acc_V = NA_real_, acc_V_hamming = NA_real_)
  jgl_precision <- c(acc_Omega = NA_real_, acc_Omega_F = NA_real_)
  jgl_lambda1 <- NA_real_
  jgl_lambda2 <- NA_real_

  jgl_fit <- tryCatch(
    fit_estimate_group_network(
      Sigma_hat_list = Sigma_hat_list,
      nvec = nvec,
      fixed_lambda1 = jgl_lambda1_fixed,
      nlambda2 = jgl_nlambda2,
      lambda2_min_ratio = jgl_lambda2_min_ratio,
      logseql2 = jgl_logseql2,
      criterion = jgl_criterion,
      gamma = jgl_gamma,
      dec = jgl_dec,
      penalty = jgl_penalty,
      weights = jgl_weights,
      maxiter = jgl_maxiter,
      truncate = jgl_truncate,
      support_tol = jgl_support_tol
    ),
    error = function(e) {
      jgl_status <<- 0
      NULL
    }
  )

  if (!is.null(jgl_fit)) {
    jgl_block <- block_metrics(V, jgl_fit$V_hat, L)
    jgl_precision <- precision_metrics(
      jgl_fit$Omega_hat,
      Omega_list,
      est_zero_tol = jgl_support_tol
    )
    jgl_lambda1 <- jgl_fit$lambda1
    jgl_lambda2 <- jgl_fit$lambda2
  }

  time_jgl <- proc.time()[["elapsed"]] - t0_jgl
  time_total <- proc.time()[["elapsed"]] - t_total

  if (show_progress) {
    cat(sprintf(
      "Finished replication %d/%d (JGL status=%d, elapsed=%.1fs)\n",
      n_rep, N_rep, jgl_status, time_total
    ))
    flush.console()
  }

  c(
    rep = n_rep,
    acc_Z = rep_acc_Z,
    ours_acc_L = ours_block[["acc_L"]],
    ours_L_hat = ours_block[["L_hat"]],
    ours_acc_V = ours_block[["acc_V"]],
    ours_acc_V_hamming = ours_block[["acc_V_hamming"]],
    ours_acc_Omega = ours_precision[["acc_Omega"]],
    ours_acc_Omega_F = ours_precision[["acc_Omega_F"]],
    jgl_status = jgl_status,
    jgl_acc_L = jgl_block[["acc_L"]],
    jgl_L_hat = jgl_block[["L_hat"]],
    jgl_acc_V = jgl_block[["acc_V"]],
    jgl_acc_V_hamming = jgl_block[["acc_V_hamming"]],
    jgl_acc_Omega = jgl_precision[["acc_Omega"]],
    jgl_acc_Omega_F = jgl_precision[["acc_Omega_F"]],
    jgl_lambda1 = jgl_lambda1,
    jgl_lambda2 = jgl_lambda2,
    time_step1 = time_step1,
    time_step2 = time_step2,
    ours_time_block = time_ours_block,
    ours_time_precision = time_ours_precision,
    jgl_time_total = time_jgl,
    time_total = time_total
  )
}

summarize_method <- function(results_df, method, prefix) {
  data.frame(
    N = N,
    J = J,
    K = K,
    L = L,
    s = s,
    N_rep = N_rep,
    method = method,
    acc_Z = mean(results_df$acc_Z),
    acc_L = mean(results_df[[paste0(prefix, "_acc_L")]], na.rm = TRUE),
    L_hat_bias = mean(results_df[[paste0(prefix, "_L_hat")]], na.rm = TRUE) - L,
    acc_V = mean(results_df[[paste0(prefix, "_acc_V")]], na.rm = TRUE),
    block_error_1_minus_ARI = 1 - mean(results_df[[paste0(prefix, "_acc_V")]], na.rm = TRUE),
    block_match_accuracy = mean(results_df[[paste0(prefix, "_acc_V_hamming")]], na.rm = TRUE),
    acc_Omega = mean(results_df[[paste0(prefix, "_acc_Omega")]], na.rm = TRUE),
    acc_Omega_F = mean(results_df[[paste0(prefix, "_acc_Omega_F")]], na.rm = TRUE),
    time_step1 = mean(results_df$time_step1),
    time_step2 = mean(results_df$time_step2),
    time_block = if (prefix == "ours") {
      mean(results_df$ours_time_block)
    } else {
      mean(results_df$jgl_time_total, na.rm = TRUE)
    },
    time_precision = if (prefix == "ours") {
      mean(results_df$ours_time_precision)
    } else {
      mean(results_df$jgl_time_total, na.rm = TRUE)
    },
    total_time_sec = if (prefix == "ours") {
      mean(results_df$time_step1 + results_df$time_step2 +
             results_df$ours_time_block + results_df$ours_time_precision)
    } else {
      mean(results_df$time_step1 + results_df$time_step2 +
             results_df$jgl_time_total, na.rm = TRUE)
    },
    opt_lambda_ind = if (prefix == "ours") opt_lambda_ind else NA_integer_,
    opt_lambda = if (prefix == "ours") opt_lambda else NA_real_,
    jgl_nlambda1 = if (prefix == "jgl") jgl_nlambda1 else NA_integer_,
    jgl_nlambda2 = if (prefix == "jgl") jgl_nlambda2 else NA_integer_,
    jgl_fixed_lambda1 = if (prefix == "jgl") jgl_lambda1_fixed else NA_real_,
    jgl_mean_selected_lambda2 = if (prefix == "jgl") mean(results_df$jgl_lambda2, na.rm = TRUE) else NA_real_,
    jgl_criterion = if (prefix == "jgl") jgl_criterion else NA_character_,
    jgl_penalty = if (prefix == "jgl") jgl_penalty else NA_character_,
    jgl_strategy = if (prefix == "jgl") jgl_strategy else NA_character_,
    jgl_status_rate = if (prefix == "jgl") mean(results_df$jgl_status) else NA_real_
  )
}

set.seed(123)

# ==============================================================================
# Define Model Parameters
# ==============================================================================
s <- env_num("SIM_S", 0.2)             # sparsity level; 0.2 or 0.6
N <- env_int("SIM_N", 1000)            # sample size
J <- env_int("SIM_J", 100)             # total number of variables
K <- env_int("SIM_K", 5)              # number of latent classes
L <- env_int("SIM_L", 5)               # number of true blocks
N_rep <- env_int("SIM_N_REP", 10)
base_seed <- env_int("SIM_SEED", 123)
nstart <- env_int("SIM_KMEANS_NSTART", 50)

M <- J / L
block_sizes <- rep(M, L)
V <- rep(seq_along(block_sizes), times = block_sizes)
pi_vec <- rep(1 / K, K)
C_vec <- c(rep(2, J / 2), rep(4, J / 2))
C_mu <- 1.0

opt_lambda_ind <- env_int("SIM_LAMBDA_IND", 2)
lambda_vec <- c(2, 1, 0.5, 0.1) * sqrt(log(J * K) / N)
if (opt_lambda_ind < 1 || opt_lambda_ind > length(lambda_vec)) {
  stop("SIM_LAMBDA_IND must be between 1 and ", length(lambda_vec), call. = FALSE)
}
opt_lambda <- lambda_vec[opt_lambda_ind]

# JGL / EstimateGroupNetwork controls. Lambda1 is fixed and lambda2 is selected
jgl_lambda1_fixed <- 0.02
jgl_nlambda1 <- NA_integer_
jgl_nlambda2 <- env_int("JGL_NLAMBDA2", 5)
jgl_lambda2_min_ratio <- env_num("JGL_LAMBDA2_MIN_RATIO", 0.01)
jgl_logseql2 <- as_bool(Sys.getenv("JGL_LOGSEQL2", unset = "true"))
jgl_criterion <- Sys.getenv("JGL_CRITERION", unset = "bic")
jgl_gamma <- env_num("JGL_GAMMA", 0.5)
jgl_dec <- env_int("JGL_DEC", 5)
jgl_strategy <- "fixed.lambda1"
jgl_penalty <- Sys.getenv("JGL_PENALTY", unset = "fused")
jgl_weights <- Sys.getenv("JGL_WEIGHTS", unset = "sample.size")
jgl_maxiter <- env_int("JGL_MAXITER", 100)
jgl_truncate <- env_num("JGL_TRUNCATE", 1e-5)
jgl_support_tol <- env_num("JGL_SUPPORT_TOL", 1e-6)
show_progress <- as_bool(Sys.getenv("SIM_PROGRESS", unset = "true"))

default_cores <- max(1, parallel::detectCores() - 1)
num_cores <- min(env_int("SIM_CORES", default_cores), N_rep)

cat(sprintf(
  "Running comparison: N=%d, J=%d, K=%d, L=%d, s=%.1f, reps=%d, cores=%d\n",
  N, J, K, L, s, N_rep, num_cores
))
cat(sprintf(
  "LCM_glasso rho=%.4g; JGL penalty=%s, criterion=%s, fixed lambda1=%.3g, nlambda2=%d\n",
  opt_lambda, jgl_penalty, jgl_criterion, jgl_lambda1_fixed, jgl_nlambda2
))

# ==============================================================================
# Replication Loop
# ==============================================================================
if (num_cores > 1) {
  cl <- parallel::makeCluster(num_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, varlist = "local_libs", envir = environment())
  parallel::clusterEvalQ(cl, {
    if (length(local_libs) > 0) {
      .libPaths(unique(c(local_libs, .libPaths())))
    }
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = setdiff(ls(envir = environment()), c("cl")),
    envir = environment()
  )
  results_list <- parallel::parLapply(cl, seq_len(N_rep), run_replication)
} else {
  results_list <- lapply(seq_len(N_rep), run_replication)
}

results_matrix <- do.call(rbind, results_list)
results_df <- as.data.frame(results_matrix)
results_df$N <- N
results_df$J <- J
results_df$K <- K
results_df$L <- L
results_df$s <- s

# ==============================================================================
# Summary Results
# ==============================================================================
summary_df <- rbind(
  summarize_method(results_df, "LCM_glasso", "ours"),
  summarize_method(results_df, "JGL_EstimateGroupNetwork", "jgl")
)

cat("\n=== Final Averaged Results ===\n")
print(summary_df[, c(
  "method",
  "acc_Z",
  "acc_L",
  "L_hat_bias",
  "block_error_1_minus_ARI",
  "block_match_accuracy",
  "acc_Omega",
  "acc_Omega_F",
  "total_time_sec"
)], row.names = FALSE)
