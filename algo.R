library(conflicted)
library(tidyverse)
library(here)
library(mclust)
library(dbscan)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

PROJ_ROOT <- here()
source(paste0(PROJ_ROOT, "/validate.R"))

# Run K-means on dataset specified by `ds`.
# Number of centers ranges from 2 to `max_num_cl`.
# For each number of centers, run K-means 10 times and select the clustering with lowest WCSS.
# Return the clustering with the highest OASW.
run_kmeans <- function(ds, max_num_cl) {
  
  dissim_matr<- dist(ds)
  
  res <- tibble(Num_cl = 2:max_num_cl) %>%
    mutate(Kmeans = map(Num_cl, ~ kmeans(ds, .x, iter.max = 1e5, nstart = 10))) %>%
    mutate(Cl = map(Kmeans, ~ .x$cluster)) %>%
    mutate(Sils = map(Cl, ~ get_sils(.x, dissim_matr))) %>%
    mutate(OASW = map_dbl(Sils, ~ .x$Sil %>% mean))
    
  opt_res <- res %>% slice_max(OASW)
  
  ret <- vector(mode = "list", length = 4)
  names(ret) <- c("num_cl", "cl", "sils", "oasw")
  ret$num_cl <- opt_res$Num_cl[1]
  ret$cl <- opt_res$Cl[[1]]
  ret$sils <- opt_res$Sils[[1]]
  ret$oasw <- opt_res$OASW[1]
  
  return(ret)
}

# Run GMMs on dataset specified by `ds`.
# Number of components ranges from 2 to `max_num_cl`.
# Try out all covariance parameterizations included in `mclust` package.
# Initialize each GMM with horizontal cut obtained from model-based hierarchical clustering.
# Return the model with the highest BIC score and the model with the highest ICL score.
run_gmms <- function(ds, max_num_cl) {
  
  # model-based hierarchical clustering to initialize components
  hc1 <- hc(ds, modelName = "VVV", use = "VARS")
  
  # Build multiple GMMs. Select the optimal model according to BIC.
  mod <- Mclust(
    ds, 
    G = 2:max_num_cl, 
    modelNames = mclust.options("emModelNames"),
    initialization = list(hcpairs = hc1),
    verbose = F
  )
  
  BIC <- NULL
  BIC$cov <- mod$modelName
  BIC$num_cl <- mod$G
  BIC$cl <- mod$classification %>% as.integer
  BIC$bic <- mod$bic
  BIC$icl <- mod$icl
  BIC$means <- array_branch(mod$parameters$mean, margin = 2)
  BIC$sigmas <- array_branch(mod$parameters$variance$sigma, margin = 3)
  
  return(BIC)
  # # Fit a GMM for each covariance parameterization and number of cluster from 2 to `max_num_cl`.
  # mods <- expand_grid(Cov = mclust.options("emModelNames"), Num_cl = 2:max_num_cl) %>%
  #   mutate(GMM = map2(Cov, Num_cl, ~ {
  #     Mclust(ds, G = .y, modelNames = .x, verbose = F, initialization = list(hcpairs = hc1))})) %>%
  #   drop_na()
  # mods <- mods %>%
  #   mutate(BIC = map_dbl(GMM, ~ .x$BIC), ICL = map_dbl(GMM, ~ .x$icl))
  # 
  # # Choose the model with the highest BIC value and the model with the highest ICL value.
  # # This model may be the 
  # mod_bic <- mods %>% slice_max(BIC, with_ties = F)
  # mod_icl <- mods %>% slice_max(ICL, with_ties = F)
  # 
  # # Create a named list to return.
  # BIC <- NULL
  # BIC$cov <- mod_bic$Cov
  # BIC$num_cl <- mod_bic$Num_cl
  # BIC$cl <- mod_bic$GMM[[1]]$classification %>% as.integer
  # BIC$bic <- mod_bic$BIC
  # BIC$icl <- mod_bic$ICL
  # BIC$means <- array_branch(mod_bic$GMM[[1]]$parameters$mean, margin = 2)
  # BIC$sigmas <- array_branch(mod_bic$GMM[[1]]$parameters$variance$sigma, margin = 3)
  # 
  # ICL <- NULL
  # ICL$cov <- mod_icl$Cov
  # ICL$num_cl <- mod_icl$Num_cl
  # ICL$cl <- mod_icl$GMM[[1]]$classification %>% as.integer
  # ICL$bic <- mod_icl$BIC
  # ICL$icl <- mod_icl$ICL
  # ICL$means <- array_branch(mod_icl$GMM[[1]]$parameters$mean, margin = 2)
  # ICL$sigmas <- array_branch(mod_icl$GMM[[1]]$parameters$variance$sigma, margin = 3)
  # 
  # ret <- NULL
  # ret$BIC <- BIC
  # ret$ICL <- ICL
  # 
  # return(ret)
}

run_hdbscan <- function(dist_obj, d, min_pts) {
  
  dissim_matr <- as.matrix(dist_obj)
  num_objs <- nrow(dissim_matr)
  
  # Single linkage hierarchy
  single_hier <-  hclust(dist_obj, method = "single")
  
  # Cuts from single linkage hierarchy.
  # Compute DBCV for each cut.
  single_cuts <- tibble(k = 2:(num_objs - 1)) %>%
    mutate(cut = map(k, ~ cutree(single_hier, k = .x))) %>%
    mutate(num_cl = map_int(cut, ~ n_distinct(.x[.x != 0]))) %>%
    filter(num_cl > 1) %>%
    mutate(val_index = map(cut, ~ get_dbcv(dissim_matr, d, .x))) %>%
    mutate(dbcv = map_dbl(val_index, ~ {
      if (is.null(.x)) {
        NA
      } else {
        sum(.x$val * .x$size_frac) 
      }
    }))
  
  # Select cut with highest DBCV. 
  opt1 <- single_cuts %>% slice_max(dbcv, with_ties = F)
  
  # Process optimal cut.
  # Handle the case where DBCV could not be computed for any cut.
  single_cut <- NULL
  if(is.na(opt1$dbcv)) {
    single_cut$min_pts <- NA
    single_cut$n <- NA
    single_cut$num_noise <- NA
    single_cut$cl <- NA
    single_cut$num_cl <- NA
    single_cut$val_index <- NA
    single_cut$dbcv <- NA
  } else {
    single_cut$min_pts <- 1
    single_cut$n <- num_objs
    opt_cut <- opt1$cut[[1]]
    single_cut$num_noise <- sum(opt_cut == 0)
    single_cut$cl <- opt_cut
    single_cut$num_cl <- opt_cut[opt_cut != 0] %>% n_distinct
    single_cut$val_index <- opt1$val_index[[1]]
    single_cut$dbcv <- opt1$dbcv
  }
  
  # # Apply FOSC to single linkage hierarchy.
  # # Compute DBCV for each partition.
  # single_foscs <- tibble(min_pts = as.integer(min_pts)) %>%
  #   mutate(cl = map(min_pts, ~ extractFOSC(single_hier, minPts = .x)$cluster)) %>%
  #   mutate(val_index = map(cl, ~ get_dbcv(dissim_matr, d, .x))) %>%
  #   mutate(dbcv = map_dbl(val_index, ~ {
  #     if (is.null(.x)) {
  #       NA
  #     } else {
  #       sum(.x$val * .x$size_frac) 
  #     }
  #   }))
  # 
  # # Select partition with highest DBCV value.
  # opt2 <- single_foscs %>% slice_max(dbcv, with_ties = F)
  # 
  # # Process optimal partition.
  # # Handle the case where DBCV could not be computed for any partition.
  # single_fosc <- NULL
  # if(is.na(opt2$dbcv)) {
  #   single_fosc$min_pts <- NA
  #   single_fosc$n <- NA
  #   single_fosc$num_noise <- NA
  #   single_fosc$cl <- NA
  #   single_fosc$num_cl <- NA
  #   single_fosc$val_index <- NA
  #   single_fosc$dbcv <- NA
  # } else {
  #   opt_cl <- opt2$cl[[1]]
  #   single_fosc$min_pts <- opt2$min_pts
  #   single_fosc$n <- num_objs
  #   single_fosc$num_noise <- sum(opt_cl == 0)
  #   single_fosc$cl <- opt_cl
  #   single_fosc$num_cl <- opt_cl[opt_cl != 0] %>% n_distinct
  #   single_fosc$val_index <- opt2$val_index[[1]]
  #   single_fosc$dbcv <- opt2$dbcv
  # }
  
  # HDBSCAN models.
  hdbscan_mods <- map(min_pts, ~ hdbscan(dist_obj, .x))
  
  # Cuts for each HDBSCAN model
  hdbscan_cuts <- NULL
  for (mod in hdbscan_mods) {
    eps_vals <- sort(mod$hc$height, decreasing = T) + .Machine$double.eps
    cuts <- tibble(min_pts = mod$minPts, eps = eps_vals) %>%
      mutate(cut = map(eps, ~ cut_tree(mod$hc, .x, kNNdist(dist_obj, mod$minPts - 1)) %>% as.integer)) %>%
      mutate(num_cl = map_int(cut, ~ n_distinct(.x[.x != 0]))) %>%
      filter(num_cl > 1)
    hdbscan_cuts <- bind_rows(hdbscan_cuts, cuts)
  }
  
  # Compute DBCV for each HDBSCAN cut.
  hdbscan_cuts <- hdbscan_cuts %>%
    mutate(val_index = map(cut, ~ get_dbcv(dissim_matr, d, .x))) %>%
    mutate(dbcv = map_dbl(val_index, ~ {
      if (is.null(.x)) {
        NA
      } else {
        sum(.x$val * .x$size_frac) 
      }
    }))
  
  # Select partition with highest DBCV value.
  opt3 <- hdbscan_cuts %>% slice_max(dbcv, with_ties = F)
  
  # Process optimal cut.
  # Handle the case where DBCV could not be computed for any cut.
  hdbscan_cut <- NULL
  if (is.na(opt3$dbcv)) {
    hdbscan_cut$min_pts <- NA
    hdbscan_cut$n <- NA
    hdbscan_cut$num_noise <- NA
    hdbscan_cut$cl <- NA
    hdbscan_cut$num_cl <- NA
    hdbscan_cut$val_index <- NA
    hdbscan_cut$dbcv <- NA
  } else {
    opt_cut <- opt3$cut[[1]]
    hdbscan_cut$min_pts <- opt3$min_pts
    hdbscan_cut$n <- num_objs
    hdbscan_cut$num_noise <- sum(opt_cut == 0)
    hdbscan_cut$cl <- opt_cut
    hdbscan_cut$num_cl <- opt_cut[opt_cut != 0] %>% n_distinct
    hdbscan_cut$val_index <- opt3$val_index[[1]]
    hdbscan_cut$dbcv <- opt3$dbcv
  }
  
  # Apply FOSC to each HDBSCAN hierarchy
  hdbscan_foscs <- tibble(
    min_pts = map_int(hdbscan_mods, ~ .x$minPts %>% as.integer),
    cl = map(hdbscan_mods, ~ .x$cluster %>% as.integer)
  ) %>%
    mutate(val_index = map(cl, ~ get_dbcv(dissim_matr, d, .x))) %>%
    mutate(dbcv = map_dbl(val_index, ~ {
      if (is.null(.x)) {
        NA
      } else {
        sum(.x$val * .x$size_frac) 
      }
    }))
  
  # Select partition with highest DBCV value.
  opt4 <- hdbscan_foscs %>% slice_max(dbcv, with_ties = F)
  
  # Process optimal partition.
  # Handle the case where DBCV could not be computed for any partition.
  hdbscan_fosc <- NULL
  if(is.na(opt4$dbcv)) {
    hdbscan_fosc$min_pts <- NA
    hdbscan_fosc$n <- NA
    hdbscan_fosc$num_noise <- NA
    hdbscan_fosc$cl <- NA
    hdbscan_fosc$num_cl <- NA
    hdbscan_fosc$val_index <- NA
    hdbscan_fosc$dbcv <- NA
  } else {
    opt_cl <- opt4$cl[[1]]
    hdbscan_fosc$min_pts <- opt4$min_pts
    hdbscan_fosc$n <- num_objs
    hdbscan_fosc$num_noise <- sum(opt_cl == 0)
    hdbscan_fosc$cl <- opt_cl
    hdbscan_fosc$num_cl <- opt_cl[opt_cl != 0] %>% n_distinct
    hdbscan_fosc$val_index <- opt4$val_index[[1]]
    hdbscan_fosc$dbcv <- opt4$dbcv
  }
  
  # return value
  ret <- NULL
  ret$single_cut <- single_cut
  ret$hdbscan_cut <- hdbscan_cut
  if(is.na(single_cut$dbcv)) {
    ret$cut <- hdbscan_cut 
  } else if(is.na(hdbscan_cut$dbcv)) {
    ret$cut <- single_cut
  } else {
    if (single_cut$dbcv >= hdbscan_cut$dbcv) {
      ret$cut <- single_cut
    } else {
      ret$cut <- hdbscan_cut
    }
  }
  ret$fosc <- hdbscan_fosc
  
  return(ret)
}

run_single_linkage <- function(dist_obj, d, min_pts) {
  
  h <-  hclust(dist_obj, method = "single")
  dissim_matr <- dist_obj %>% as.matrix
  
  h_cl <- extractFOSC(h, minPts = min_pts)$cluster
  h_val_index <- get_dbcv(dissim_matr, d, h_cl)
  if (is.null(h_val_index)) {
    h_val_index <- NA
    h_dbcv <- NA
  } else {
    h_dbcv <- sum( h_val_index$val * h_val_index$size_frac )
  }
  
  fosc_cut <- NULL
  fosc_cut$min_pts <- min_pts
  fosc_cut$n <- nrow(dissim_matr)
  fosc_cut$cl <- h_cl
  fosc_cut$num_cl <- h_cl[h_cl != 0] %>% n_distinct
  fosc_cut$num_noise <- sum(h_cl == 0)
  fosc_cut$val_index <- h_val_index
  fosc_cut$dbcv <- h_dbcv
  
  max_num_cl <- nrow(dissim_matr) - 1
  horiz_cuts <- tibble(Num_cl = 2:max_num_cl) %>%
    mutate(Cut = map(Num_cl, ~ cutree(h, k = .x))) %>%
    mutate(Num_cl = map_int(Cut, ~ n_distinct(.x[.x != 0]))) %>%
    filter(Num_cl > 1)
  
  horiz_cuts <- horiz_cuts %>%
    mutate(Val_index = map(Cut, ~ get_dbcv(dissim_matr, d, .x))) %>%
    mutate(DBCV = map_dbl(Val_index, ~ {
      if (is.null(.x)) {
        NA
      } else {
        sum(.x$val * .x$size_frac)
      }
    }))
  
  opt_cut <- NULL
  opt_cut$min_pts <- min_pts
  opt_cut$n <- nrow(dissim_matr)
  # If all horizontal cuts cannot be evaluated by DBCV
  if (all(map_lgl(horiz_cuts$DBCV, is.na))) {
    opt_cut$cl <- NA
    opt_cut$num_cl <- NA
    opt_cut$num_noise <- NA
    opt_cut$val_index <- NA
    opt_cut$dbcv <- NA
  } else {
    opt_horiz_cut <- horiz_cuts %>% slice_max(DBCV, with_ties = F)
    opt_cl <- opt_horiz_cut$Cut[[1]]
    opt_cut$cl <- opt_cl
    opt_cut$num_cl <- opt_horiz_cut$Num_cl
    opt_cut$num_noise <- sum(opt_cl == 0)
    opt_cut$val_index <- opt_horiz_cut$Val_index[[1]]
    opt_cut$dbcv <- opt_horiz_cut$DBCV
  }
  
  ret <- NULL
  ret$FOSC <- fosc_cut
  ret$HORIZ <- opt_cut
  
  return(ret)
}

# cutree from `stats` package does not distinguish noise as 0, so we make a new
# method to do it manually.
# See https://cran.r-project.org/web/packages/dbscan/vignettes/hdbscan.html
cut_tree <- function(hcl, eps, core_dist) {
  
  cuts <- unname(cutree(hcl,h = eps))
  cuts[which(core_dist > eps)] <- 0 # Use core distance to distinguish noise
  
  return(cuts)
}
