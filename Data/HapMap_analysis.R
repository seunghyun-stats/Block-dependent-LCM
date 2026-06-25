library(ggplot2)
library(GGally)
library(reshape2)
library(RSpectra)
library(mvtnorm)
library(RcppHungarian)
library(polycor)
library(dplyr)
library(igraph)
library(huge)

perm <- function(x, p, K) {
  x_perm <- rep(NA_integer_, length(x))
  for (k in seq_len(K)) {
    x_perm[x == k] <- p[k]
  }
  x_perm
}

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

#### prepare data ####
setwd("C:/Users/Sky/Documents/HapMap")
path_root <- "C:/Users/Sky/Documents/Hapmap/Archive/" 
# raw data available at https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3

# obtain genetic and ethnicity data
data_hapmap <- read.csv(paste0(path_root, "merged_add.traw"), sep="\t")
data_gene <- t(data_hapmap[, -(1:6)])
data_hapmap <- data_hapmap[, 1:6]

eth <- read.csv(paste0(path_root, "hapmap_eth.csv"))[, 2]
eth_num <- as.numeric(as.factor(eth))

# consider a subset of 4 ethnicity groups
subset <- which(eth %in% c('CEU', 'YRI', 'MEX', 'GIH'))
cluster_true = eth[subset]
cluster_true = as.numeric(as.factor(cluster_true))


# select the location with larger variance
data_gene_22 = data_gene[subset, which(data_hapmap[, 1] == 22)]
var_vec = apply(data_gene_22,2,var)
top_500_indices <- sort(order(var_vec, decreasing = TRUE)[1:500])
R_mat = data_gene_22[, top_500_indices]

data_hapmap_updated = data_hapmap[which(data_hapmap[, 1] == 22), ]
data_hapmap_500 = data_hapmap_updated[top_500_indices, ] # hapmap file corresponding to selected
data_hapmap_500[1:176,]


# data matrix R (takes values 0, 1, 2.)
N <- nrow(R_mat)
J <- ncol(R_mat)
K <- 4

###### Step 1: Clustering ######
svd_res <- svds(R_mat, K)

kmeans_res <- kmeans(svd_res$u, centers=K, iter.max=100, nstart=10)
Z_hat <- kmeans_res$cluster
centers <- kmeans_res$centers

# find the optimal permutation
mat <- matrix(0, K, K)
for (i in 1:K) {
  for (j in 1:K) {
    mat[i, j] <- sum((Z_hat == i) & (cluster_true == j))
  }
}

perm_mat <- HungarianSolver(-mat)$pairs[, 2]
Z_hat <- perm(Z_hat, perm_mat, K)
sum(Z_hat != cluster_true) # Hamming clustering error


###### Step 2: Covariance matrix estimation ###### 
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
  
  pc_cor <- hetcor(R_k, ML = FALSE, std.err = FALSE)$correlations
  Sigma_k <- pc_cor
  Sigma_k[is.na(Sigma_k)] <- 0
  Sigma_hat_list[[k]] <- Sigma_k
}

##### Block structure recovery #####
sum_I <- matrix(0, nrow = J, ncol = J)
for (k in 1:K) {
  sum_I <- sum_I + abs(Sigma_hat_list[[k]]) 
}
diag(sum_I) <- 0
avg_support <- sum_I / K

g <- graph_from_adjacency_matrix(avg_support, mode = "undirected", weighted = TRUE)

# default resolution
res = 1
set.seed(2026)
leiden_res <- cluster_leiden(g, "modularity", resolution = res)
L_hat_low <- length(leiden_res)
V_hat_low <- membership(leiden_res)

# higher resolution
res = 1.5
leiden_res <- cluster_louvain(g,resolution = res)
L_hat_high <- length(leiden_res)
V_hat_high <- membership(leiden_res)

##### compute the average off-block correlation #####
avg_abs_offblock <- function(Sigma_list, V_hat) {
  J <- length(V_hat)
  
  # unordered off-diagonal pairs not in the same block
  offblock_mask <- outer(V_hat, V_hat, "!=") & upper.tri(matrix(TRUE, J, J))
  
  mean(unlist(lapply(Sigma_list, function(Sigma_k) {
    abs(Sigma_k[offblock_mask])
  })), na.rm = TRUE)
}

avg_offblock_low  <- avg_abs_offblock(Sigma_hat_list, V_hat_low)
avg_offblock_high <- avg_abs_offblock(Sigma_hat_list, V_hat_high)

cat("Average |sigma| outside low-resolution blocks:",
    avg_offblock_low, "\n")

cat("Average |sigma| outside high-resolution blocks:",
    avg_offblock_high, "\n")
