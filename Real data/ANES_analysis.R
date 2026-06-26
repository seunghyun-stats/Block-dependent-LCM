library(reshape2)
library(ggplot2)
library(mclust)
library(MASS)
library(Matrix)
library(tictoc)
library(RSpectra)
library(RcppHungarian)
library(igraph)
library(polycor)       
library(huge)          
library(qgraph)

perm <- function(x, p, K) {
  x_perm <- rep(NA, length(x))
  for (i in 1:length(x)) {
    for (k in 1:K) {
      x_perm[x == k] <- p[k]
    }
  }
  return(x_perm)
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

setwd("C:/Users/Sky/Documents")
R_sub = read.csv("R_anes.csv")[, -1]
# the raw dataset is available at https://electionstudies.org/data-center/
# the processing follows the R code ANES.R available at https://github.com/lscientific/gGoM/blob/main/real_data/ANES.R


### Record the baseline block structure (available at the questionnaire at link) ###
colnames(R_sub)
block_true = list(1,2,3,4:10,11:16,17:19,20:22,23:30,31:44,45:56,57,58,59:62,63:64,c(65:71,137:138),72:77,
                  78:80,81:86,87:95,96:107,108,109,110,111:113,114:124,125:126,127:131,132,133:136,
                  139:140,141:144
)
names(block_true) = c("FOLLOW POLITICS","VOTER REGISTRATION","TURNOUT","PARTICIPATION",
                      "GLOBAL EMOTION","PRESIDENTIAL APPROVAL","ECONOMIC PERFORMANCE","INFLATION",
                      "ISSUE IMPORTANCE", "Issue ownership", "Climate", "Ukraine", "Trust Experts",
                      "Political disagreement", "ABORTION", "ABORTION EMOTIONS", "GUNS AND CRIME", 
                      "IMMIGRANT EMOTIONS", "DEMOCRATIC ATTITUDES", "ELECTORAL INTEGRITY", "POLITICAL EFFICACY",
                      "RACISM", "FEMINIST ATTITUDES", "POLITICAL TOLERANCE", "RACIAL STEREOTYPES", "IDENTITIES",
                      "ROLE OF SCHOOLS", "GREAT REPLACEMENT", "RACIAL PRIVILEGE", "TRANSGENDER ATTITUDES", "RACIAL RESENTMENT"
)
V_true <- numeric(max(unlist(block_true)))

for (i in seq_along(block_true)) {
  V_true[block_true[[i]]] <-i # names(block_true)[i] # i
}


#### estimation ####
K <- 3
R_sub = R_sub[,-137]

C_list <- apply(R_sub, 2, function(x) {
  length(unique(x))
})

N <- nrow(R_sub)
J <- ncol(R_sub)

# ==============================================================================
# Step 1: Spectral Clustering
# ==============================================================================
set.seed(2026)
R_mat = R_sub

R_flat <- flatten_ordinal(R_mat, C_list)
nstart <- 50
U <- svds(R_flat, k=K)$u
# U <- svds(as.matrix(R_mat), k=K)$u
kmeans_res <- kmeans(U, centers=K, iter.max=100, nstart=nstart)
Z_hat <- kmeans_res$cluster


### Draw Figure 4 ###
Delta_hat <- vector("list", J)
for (j in 1:J) {
  C_j <- C_list[j] # Total number of categories for item j
  
  thresh_mat <- matrix(NA, nrow = K, ncol = C_j - 1)
  rownames(thresh_mat) <- paste0("Class_", 1:K)
  colnames(thresh_mat) <- paste0("Delta_", 1:(C_j - 1))
  
  for (k in 1:K) {
    R_kj <- R_mat[Z_hat == k, j]
    N_k <- length(R_kj)
   
    counts <- table(factor(R_kj, levels = 1:C_j))
    cum_probs <- cumsum(counts) / N_k
    cum_probs <- cum_probs[-C_j]
    cum_probs <- pmax(pmin(cum_probs, 1 - 1e-5), 1e-5)
    thresh_mat[k, ] <- qnorm(cum_probs)
  }
  
  Delta_hat[[j]] <- thresh_mat
}

Delta_1_mat <- matrix(0, nrow = J, ncol = 3)
rownames(Delta_1_mat) <- colnames(R_sub)

for (j in 1:J) {
  Delta_1_mat[j, ] <- Delta_hat[[j]][, 1]
}

T_diff <- apply(Delta_1_mat, 1, function(x) {
  max(c(abs(x[1] - x[2]), abs(x[1] - x[3]), abs(x[2] - x[3])))
})

selected_indices <- order(T_diff, decreasing = TRUE)[1:25]
Delta_1_mat[selected_indices, ]
top_7_items <- Delta_1_mat[selected_indices[c(1:5,11,17)], ]

custom_labels <- c(
  "Vote for Trump in 2024? (4)", 
  "Is Biden the legitimate winner of 2020? (2)", 
  "Do Democrats better handle employment? (3)", 
  "Do Democrats better handle cost of living? (3)", 
  "Do you believe in Trump's statement \"election was stolen...\"? (2)",
  "Is ineligible voters being allowed to vote a problem? (2)", 
  "Is the Jan 6 incident a justified protest? (2)"
)

plot_data <- as.data.frame(top_7_items)
colnames(plot_data) <- c("Class_1", "Class_2", "Class_3")
plot_data$Item <- custom_labels
plot_data$Item <- factor(plot_data$Item, levels = rev(custom_labels))

plot_data_long <- melt(plot_data, id.vars = "Item", 
                       variable.name = "Class", value.name = "Delta_1")

p <- ggplot(plot_data_long, aes(x = Delta_1, y = Item, color = Class)) +
  geom_line(aes(group = Item), color = "gray80", linewidth = 1) +
  geom_point(size = 5, alpha = 0.9) + 
  theme_minimal() +
  labs(x = expression(Delta["j,1,k"]),
       y = "") +
  
  scale_color_manual(
    values = c("Class_1" = "firebrick",   # Democrat
               "Class_2" = "dodgerblue",  # Republican
               "Class_3" = "gray50"),      # Independent
    labels = c("Class_1" = "Republican (k=1)", 
               "Class_2" = "Democrat (k=2)", 
               "Class_3" = "Independent (k=3)")
  ) +
  
  theme(axis.text.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold"),
        panel.grid.major.y = element_blank())

ggsave("Top_10_Delta1_Differences.pdf", plot = p, width = 18, height = 6, device = "pdf")


# ==============================================================================
# Step 2: Covariance Matrix Estimation
# ==============================================================================
Sigma_hat_list <- list()

for (k in 1:K) {
  indices_k <- which(Z_hat == k)
  R_k <- R_mat[indices_k, , drop = FALSE]
  N_k <- length(indices_k)
  
  pc_cor <- hetcor(R_k, ML = FALSE, std.err = FALSE)$correlations
  
  Sigma_k <- pc_cor
  Sigma_k[is.na(Sigma_k)] <- 0
  Sigma_hat_list[[k]] <- Sigma_k
}

sum_I <- matrix(0, nrow = J, ncol = J)
for (k in 1:K) {
  sum_I <- sum_I + abs(Sigma_hat_list[[k]]) 
}
diag(sum_I) <- 0
avg_support <- sum_I / K

g <- graph_from_adjacency_matrix(avg_support, mode = "undirected", weighted = TRUE)

# Run Leiden community detection
res = 1
leiden_res <- cluster_leiden(g, "modularity", resolution = res)

L_hat_low <- length(leiden_res)
V_hat_low <- membership(leiden_res)
adjustedRandIndex(V_true, V_hat_low)

### higher resolution
res = 1.5
leiden_res <- cluster_leiden(g, "modularity", resolution = res)

L_hat_high <- length(leiden_res)
V_hat_high <- membership(leiden_res)
adjustedRandIndex(V_true, V_hat_high)

perm_indices <- order(V_hat_low, V_hat_high)
avg_support_permuted <- avg_support[perm_indices, perm_indices]

# Sort the labels to match the permuted matrix
V_low_sorted  <- V_hat_low[perm_indices]
V_high_sorted <- V_hat_high[perm_indices]


### Draw Figure 5 ###
bounds_low_internal  <- cumsum(rle(V_low_sorted)$lengths) + 0.5
bounds_high_internal <- cumsum(rle(V_high_sorted)$lengths) + 0.5
bounds_high_internal <- bounds_high_internal[-9] 

J <- length(perm_indices)
bounds_low  <- c(0.5, bounds_low_internal[-length(bounds_low_internal)], J + 0.5)
bounds_high <- c(0.5, bounds_high_internal[-length(bounds_high_internal)], J + 0.5)

gaps <- which(diff(bounds_high) == 1)
drop_idx <- ifelse(bounds_high[gaps + 1] %in% bounds_low, gaps, gaps + 1)
bounds_high <- if(length(drop_idx) > 0) bounds_high[-drop_idx] else bounds_high
avg_support_permuted <- avg_support[perm_indices,perm_indices]

block_edges <- c(0.5, bounds_low, J + 0.5)
block_centers <- (block_edges[-length(block_edges)] + block_edges[-1]) / 2
centers_to_label <- block_centers[-1]

z_min <- min(avg_support_permuted, na.rm = TRUE)
z_max <- max(avg_support_permuted, na.rm = TRUE)

col_palette <- colorRampPalette(c("white", "navy"))(100)

pdf("multi_resolution_comparison.pdf", width = 14, height = 7)
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(5, 5, 1))
par(mar = c(6, 5, 6, 2) + 0.1)

# Plot 1: Low Resolution
image(x = 1:J, y = 1:J, z = avg_support_permuted, col = col_palette, 
      main = expression(paste("Default Leiden (", hat(L), " = 7)")),
      xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes = FALSE,
      cex.main = 2, cex.lab = 2, cex.axis = 2)

abline(v = bounds_low, col = "black", lwd = 2)
abline(h = bounds_low, col = "black", lwd = 2)

axis(side = 1,
     at = centers_to_label,
     labels = 1:length(centers_to_label),
     tick = FALSE,
     line = 0.5,
     cex.axis = 2,
     font = 2)

# Plot 2: High Resolution
image(x = 1:J, y = 1:J, z = avg_support_permuted, col = col_palette, 
      main = expression(paste("High-resolution (", hat(L), " = 16)")),
      xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes = FALSE,
      cex.main = 2, cex.lab = 2, cex.axis = 2)

for (i in 1:(length(bounds_low) - 1)) {
  low_start <- bounds_low[i]
  low_end   <- bounds_low[i+1]
  high_in_block <- bounds_high[bounds_high > low_start & bounds_high < low_end]
  
  if (length(high_in_block) > 0) {
    high_edges <- sort(unique(c(low_start, high_in_block, low_end)))
    
    for (j in 1:(length(high_edges) - 1)) {
      rect(xleft = high_edges[j], ybottom = high_edges[j],
           xright = high_edges[j+1], ytop = high_edges[j+1],
           border = "red", lwd = 2)
    }
  }
}

abline(v = bounds_low, col = "black", lwd = 2)
abline(h = bounds_low, col = "black", lwd = 2)

axis(side = 1,
     at = centers_to_label,
     labels = 1:length(centers_to_label),
     tick = FALSE,
     line = 0.5,
     cex.axis = 2,
     font = 2)

# Plot 3: Color Bar
par(mar = c(6, 0.5, 6, 5.5) + 0.1)

color_seq <- seq(z_min, z_max, length.out = 100)
image(x = 1, y = color_seq, z = matrix(color_seq, nrow = 1), 
      col = col_palette, axes = FALSE, xlab = "", ylab = "", 
      main = "", cex.main = 2)

axis(4, las = 1, cex.axis = 2)
box()

dev.off()


# ==============================================================================
# Step 3: Precision Matrix Estimation
# ==============================================================================
V_hat <- V_hat_low 
# V_hat <- V_hat_high
L_labels <- unique(V_hat)
opt_lambda <-  0.3* sqrt(log(J*K)/N) # selected by CV

Omega_hat_final <- vector("list", length = K)
for (k in 1:K) Omega_hat_final[[k]] <- matrix(0, nrow = J, ncol = J)

for (k in 1:K) {
  Sigma_k <- Sigma_hat_list[[k]][perm_indices, perm_indices]
  
  for (l_idx in seq_along(L_labels)) {
    l <- L_labels[l_idx]
    indices <- which(V_low_sorted == l)
    
    if (length(indices) == 0) next
    
    Sigma_sub <- Sigma_k[indices, indices, drop = FALSE]
    
    # Ensure Positive Definiteness
    Sigma_sub <- (Sigma_sub + t(Sigma_sub)) / 2
    min_eig <- min(eigen(Sigma_sub, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig < 1e-4) {
      Sigma_sub <- Sigma_sub + diag(abs(min_eig) + 1e-4, nrow(Sigma_sub))
    }
    
    # graphical lasso
    out_final <- huge(Sigma_sub, method = "glasso", lambda = opt_lambda, verbose = FALSE)
    Omega_hat_final[[k]][indices, indices] <- as.matrix(out_final$icov[[1]])
  }
}

omega_array <- simplify2array(Omega_hat_final)

### Draw Figure 6 ###
V_sorted = V_low_sorted
# block = 6
# ind_block = which(V_sorted == block)
ind_block = c(121:123,126:133)

Omega_plot <- lapply(Omega_hat_final[1:3], function(mat) {
  diag(mat) <- pmin(diag(mat), 2)
  return(mat)
})
Omega_plot <- lapply(Omega_plot, function(mat) {
  mat[abs(mat) < 3*opt_lambda[1]] <- 0 # thresholding for better visualization
  return(mat)
})
col_palette_omega <- colorRampPalette(c("firebrick", "white", "navy"))(101)

# compute the list of partial correlations (rho)
adj_list <- lapply(Omega_plot[1:3], function(mat) {
  sub_mat <- mat[ind_block, ind_block, drop = FALSE]
  sub_mat <- (sub_mat + t(sub_mat)) / 2
  adj <- -cov2cor(sub_mat) 
  diag(adj) <- 0  
  
  return(adj)
})


max_edge <- max(abs(unlist(adj_list)), na.rm = TRUE)
same_sign <- (sign(adj_list[[1]]) == sign(adj_list[[2]])) & 
  (sign(adj_list[[2]]) == sign(adj_list[[3]])) &
  (adj_list[[1]] != 0)

# shared part of the three partial correlation networks
shared_adj <- matrix(0, nrow = nrow(adj_list[[1]]), ncol = ncol(adj_list[[1]]))
shared_adj[same_sign] <- (adj_list[[1]][same_sign] + 
                            adj_list[[2]][same_sign] + 
                            adj_list[[3]][same_sign]) / 3

avg_adj <- (abs(adj_list[[1]]) + abs(adj_list[[2]]) + abs(adj_list[[3]])) / 3
fixed_layout <- qgraph(avg_adj, layout = "spring", DoNotPlot = TRUE)$layout
node_labels <- paste0(1:11)

# code for plot
pdf("ANES_precision_network.pdf", width = 24, height = 6)
par(mfrow = c(1, 4))

# Plot the 3 class-specific networks
for (k in 1:3) {
  qgraph(adj_list[[k]], 
         layout = fixed_layout,
         directed = FALSE, 
         title = bquote(hat(Omega)[.(k)]),
         labels = node_labels,
         maximum = max_edge,
         
         posCol = "navy", 
         negCol = "firebrick",   
         fade = FALSE,  
         
         vsize = 10,
         borders = TRUE,
         title.cex = 3, 
         label.cex = 2.2)
}

# Plot the shared structure
qgraph(shared_adj, 
       layout = fixed_layout,         
       directed = FALSE,              
       title = "Shared",
       labels = node_labels,          
       maximum = max_edge,            
       posCol = "navy",               
       negCol = "firebrick",                
       fade = FALSE,                   
       vsize = 10,                    
       borders = TRUE,                
       title.cex = 3,               
       label.cex = 2.2)

dev.off()
par(mfrow = c(1, 1))

