library(conflicted)
library(tidyverse)
library(cluster)
library(mclust)
library(igraph)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("crossing", "tidyr")
conflict_prefer("as_data_frame", "igraph")

# Function to compute the silhouette values based on a dissimilarity matrix and a set of labels
get_sils <- function(lbls, dist_obj) {
  
  sils <- silhouette(lbls, dist_obj)
  
  sils1 <- tibble(
    Cluster = sils[, 1] %>% as.integer,
    Neighbor = sils[, 2] %>% as.integer,
    Sil = sils[, 3]
  )
  
  return(sils1)
}

# Function to compute silhouette values based on labels, a dataset, and a covariance matrix
# for each label.
# Note that the labels should be integers: 1,2,3, ..., m.
# The covariance matrices should be passed as a list of matrices.
get_sils_mahal <- function(lbls, ds, sigmas) {
  
  ds <- ds %>% as.matrix
  
  # Compute Mahalanobis distance between every pair of points.
  dists <- expand_grid(x = 1:nrow(ds), y = 1:nrow(ds)) %>%
    filter(x != y) %>%
    mutate(Cluster = map_int(x, ~ lbls[.x] %>% as.integer)) %>%
    mutate(Neighbor = map_int(y, ~ lbls[.x] %>% as.integer)) %>%
    mutate(Mahal = map2_dbl(x, y, ~ mahalanobis(ds[.x,], ds[.y,], cov = sigmas[[lbls[.y]]]) %>% sqrt))
  
  # For each point, compute the mean distance between the point and every cluster.
  mean_dists <- dists %>%
    group_by(x, Cluster, Neighbor) %>%
    summarize(Mean_dist = mean(Mahal), .groups = "drop")
  
  # Function to compute silhouette value
  sil <- function(df) {
    
    # compute a
    a_val <- df %>%
      filter(Cluster == Neighbor) %>%
      select(a = Mean_dist)
    
    # compute b
    b_val <- df %>%
      filter(Cluster != Neighbor) %>%
      arrange(Mean_dist) %>%
      slice_head() %>%
      select(everything(), b = Mean_dist)
    
    # Compute the point's silhouette. If the point belongs to a singleton cluster, then `a_val` will
    # have zero rows. In this case, the point's silhouette is defined to be 0.
    if (nrow(a_val) == 0) {
      sil_val = b_val %>%
        mutate(Sil = 0) %>%
        select(-b)
    } else {
      sil_val <- bind_cols(b_val, a_val) %>%
        mutate(Sil = map2_dbl(a, b, ~ (.y - .x) / max(.x, .y))) %>%
        select(-a, -b)
    }
    
    return(sil_val)
  }
  
  # Compute each point's silhouette.
  sils <- mean_dists %>%
    nest(data = c(Cluster, Neighbor, Mean_dist)) %>%
    mutate(sil = map(data, sil)) %>%
    select(-data) %>%
    unnest(sil)
  
  return(sils)
}

# Compute DBCV based on a dissimilarity matris, the dimensionality of the dataset and labels.

get_dbcv <- function(dissim_matr, d, lbls) {
  
  # DBCV operates on non-noise objects. So we delete noise objects but
  # keep track of the total number of objects (noise and non-noise).
  num_objs <- nrow(dissim_matr)
  dissim_matr <- dissim_matr[lbls != 0, lbls != 0]
  lbls <- lbls[lbls != 0]
  
  if (length(lbls) == 0) {
    # print("Every point is noise!")
    return(NULL)
  }
  
  # If there is a singleton cluster, the object's all-points-core-distance
  # is undefined. Therefore, return NULL
  if (1 %in% table(lbls)) {
    # print("Singleton cluster!")
    return(NULL)
  }
  
  # If there is only one set within the partition, then cannot compute
  # density separation. Therefore, return NULL
  if (n_distinct(lbls) == 1) {
    return(NULL)
  }
  
  # compute core distance of each object
  core_dists <- map_dbl(1:nrow(dissim_matr), ~ all_points_core_dist(.x, dissim_matr, d, lbls))
  
  # compute mutual reachability distance matrix.
  mrd <- mutual_reachability_distance(dissim_matr, core_dists)
  rownames(mrd) <- 1:nrow(dissim_matr)
  colnames(mrd) <- 1:nrow(dissim_matr)
  
  # compute mrd graph
  mrd_graph <- graph_from_adjacency_matrix(mrd, weighted = T, mode = "upper")
  
  # compute density sparseness of each cluster
  lbls_unique <- lbls %>% unique %>% sort
  DSC <- tibble(cl = lbls_unique) %>%
    mutate(val = map_dbl(cl, ~ density_sparseness(mrd_graph, .x, lbls)))
  
  if (any(map_lgl(DSC$val, is.na))) {
    return(NULL)
  }
  
  # compute density separation for each pair of clusters
  label_combs <- combn(lbls_unique, 2, simplify = F)
  DSPC <- tibble(
    cl1 = map_int(label_combs, ~ .x[1] %>% as.integer),
    cl2 = map_int(label_combs, ~ .x[2] %>% as.integer)
  ) %>%
    mutate(val = map2_dbl(cl1, cl2, ~ density_separation(mrd, mrd_graph, .x, .y, lbls)))
  
  if (any(map_lgl(DSPC$val, is.na))) {
    return(NULL)
  }
  
  # compute validity index and fractional size for each cluster
  val_index <- tibble(cl = lbls_unique) %>%
    mutate(val = map_dbl(cl, ~ validity_index_cluster(.x, DSPC, DSC))) %>%
    mutate(size_frac = map_dbl(cl, ~  sum(lbls == .x)/ num_objs))
  
  return(val_index)
}

# compute all-points-core-dist for an object
all_points_core_dist <- function(obj, dissim_matr, d, labels) {
  
  # get distances to every object in the same cluster as i. 
  dists <- dissim_matr[obj, labels == labels[obj]]
  
  # compute numerator
  numer_terms <- (1 / dists[dists != 0]) ^ d
  numer <- sum(numer_terms)
  
  core_dists <- (numer / (length(dists) - 1)) ^ (-1 / d)
  return(core_dists)
}

# compute mutual reachability distance matrix
mutual_reachability_distance <- function(dissim_matr, core_dists) {
  
  num_objects <- nrow(dissim_matr)
  mrd <- matrix(rep(0, num_objects^2), nrow = num_objects)
  for (i in 1:num_objects) {
    
    j <- i + 1
    while(j <= num_objects) {
      mrd[i, j] <- max(dissim_matr[i,j], core_dists[i], core_dists[j])
      j <- j + 1
    }
  }
  return(mrd)
}

# compute density sparseness
density_sparseness <- function(mrd_graph, cl, labels) {
  
  mst_cl <- induced_subgraph(mrd_graph, which(labels == cl)) %>% mst
  degs <- degree(mst_cl)
  es <- mst_cl %>% as_data_frame
  
  internal_es <- map_lgl(1:nrow(es), ~ getElement(degs, es[.x, 1]) != 1 && getElement(degs, es[.x, 2]) != 1)
  
  if (!any(internal_es)) {
    # print(paste("Density sparseness of cluster", cl, "cannot be computed because the MST of cluster",
    #             cl, "does not have internal edges."))
    return(NA)
  } else {
    max_internal_e <- es[internal_es,]$weight %>% max
    return(max_internal_e)
  }
}

# compute density separation between two clusters
density_separation <- function(mrd, mrd_graph, cl1, cl2, labels) {
  
  mst_cl1 <- induced_subgraph(mrd_graph, which(labels == cl1)) %>% mst
  mst_cl2 <- induced_subgraph(mrd_graph, which(labels == cl2)) %>% mst
  
  degs1 <- degree(mst_cl1)
  degs2 <- degree(mst_cl2)
  
  inodes1 <- names(which(degs1 != 1)) %>% as.integer
  inodes2 <- names(which(degs2 != 1)) %>% as.integer
  
  if (length(inodes1) == 0) {
    # print(paste("Density separation between clusters", cl1, "and", cl2,
    #             "cannot be computed because the MST of", cl1, "does not have any internal nodes."))
    return(NA)
  }
  
  if (length(inodes2) == 0) {
    # print(paste("Density separation between clusters", cl1, "and", cl2,
    #             "cannot be computed because the MST of", cl2, "does not have any internal nodes."))
    return(NA)
  }
  
  dists <- expand_grid(inodes1, inodes2) %>%
    mutate(dist = map2_dbl(inodes1, inodes2, ~ mrd[min(.x, .y), max(.x, .y)]))
  
  dense_separ <- dists$dist %>% min()
  return(dense_separ)
}

# compute validity index of a cluster
validity_index_cluster <- function(cl_in, DSPC, DSC) {
  
  dsc <- DSC %>% filter(cl == cl_in) %>% pull(val) %>% unlist
  dspc <- DSPC %>% filter(cl1 == cl_in | cl2 == cl_in) %>% pull(val) %>% unlist %>% min
  
  numer <- dspc - dsc
  denom <- max(dspc, dsc)
  return(numer / denom)
}

get_ARI <- function(gt, cl) {
  return(adjustedRandIndex(gt, cl))
}

external_eval <- function(gt, cl) {
  
  # confusion matrix
  cm <- as.matrix(table(Cluster = cl, Ground_truth = gt))
  
  rsums <- apply(cm, 1, sum)
  csums <- apply(cm, 2, sum)
  
  # precision matrix 
  prec_matr <- cm / rsums
  
  # recall matrix
  rec_matr <- t( t(cm) / csums)
  
  # F-measure matrix
  f_matr <- 2 * prec_matr * rec_matr / (prec_matr + rec_matr)
  
  # maximum F-measure for each cluster
  maxfmeas <- tibble(Cluster = rownames(f_matr)) %>%
    mutate(Fmeas = map(Cluster, ~ {
      # complex(es) with the largest F-measure
      c <- names(which(f_matr[.x,] == max(f_matr[.x,], na.rm = T)));
      fmeas <- f_matr[.x, c];
      tibble(Complex = c, Fmeas = fmeas)
    })) %>%
    unnest(Fmeas) %>%
    mutate(Cluster = Cluster %>% as.integer, Complex = Complex %>% as.integer) %>%
    arrange(Cluster)
  
  # return data
  ret <- vector(mode = "list", length = 5)
  names(ret) <- c("cm", "prec", "rec", "fmeas", "maxfmeas")
  ret$cm <- cm
  ret$prec <- prec_matr
  ret$rec <- rec_matr
  ret$fmeas <- f_matr
  ret$maxfmeas <- maxfmeas
  
  return(ret)
}

# set.seed(26)
# d <- bind_rows(
#   tibble(x = runif(3, min= .25, max = .75), y = runif(3, min = .25, max = .75)),
#   tibble(x = runif(3, min = -.75, max = -.25), y = runif(3, min = .25, max = .75)),
#   tibble(x = runif(3, min= .25, max = .75), y = runif(3, min = -.75, max = -.25))
# )
# 
# ggplot(d) +
#   geom_point(aes(x = x, y = y)) +
#   expand_limits(x = c(-1, 1), y = c(-1, 1))
# 
# dissim_matr <- dist(d) %>% as.matrix
# d <- 2
# lbls <- c(rep(1L, 3), rep(2L, 3), rep(3L, 3))
