library(conflicted)
library(tidyverse)
library(here)
library(progress)
library(proxy)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")
conflict_prefer("dist", "proxy")

source(paste0(here(), "/load_data.R"))
source(paste0(here(), "/algo.R"))
source(paste0(here(), "/validate.R"))

# Load melting data from K562 intact cell line.
mdata <- load_mdata()

# Load CORUM protein complex data in long format
cdata_long <- load_cdata()

# Remove proteins from `cdata_long` with no melting data, as we cannot analyze these.
cdata_long <- cdata_long %>% semi_join(mdata, by = "Protein")

# `cdata` contains the same data as `cdata_long` but not in long format
cdata <- cdata_long %>%
  nest(Prots = Protein) %>%
  mutate(Prots = map(Prots, ~ unlist(.x, use.names = F)))
ids <- cdata$Complex
cdata <- cdata$Prots
names(cdata) <- ids

# Remove complexes with fewer than three proteins. Such complexes may be present because
# we have removed proteins with no melting data
cdata <- cdata %>% keep(~ length(.x) >= 3)
cdata_long <- cdata_long %>% filter(Complex %in% (names(cdata) %>% as.integer))

# Remove proteins from melting data that are not present in the complex data
mdata <- mdata %>% semi_join(cdata_long, by = "Protein")

# Fit a three parameter log-logistic model based on each protein's melting data.
mparams <- param_mdata(mdata)

# Function to choose `num` non-overlapping complexes. The first complex is the `seed`. 
# The other `num` - 1 complexes are chosen at random.
choose_complexes <- function(num, seed) {
  
  # choose the seed complex
  chosen_ids <- seed
  chosen_prots <- cdata[[seed]]
  
  # Remove seed complex so that it is not chosen again
  cdata[[seed]] <- NULL
  cdata <- cdata %>% compact
  
  # choose the other complexes
  for (i in 1:(num - 1)) {
    
    repeat {
      
      # choose a complex at random
      a_complex <- sample(cdata, 1)
      id <- names(a_complex)
      prots <- a_complex[[1]]
      
      # Break if the complex does not have any proteins in common with the already chosen complexes
      if (length(intersect(prots, chosen_prots)) == 0) break
    }
    
    # Add the complex to those already chosen
    chosen_ids <- c(chosen_ids, id)
    chosen_prots <- c(chosen_prots, prots)
    
    # Remove the complex so that it is not chosen again
    cdata[[seed]] <- NULL
    cdata <- cdata %>% compact
  }
  
  return(cdata_long %>% filter(Complex %in% chosen_ids))
}

# Repeatedly choose samples of non-overlapping complexes.
cdata_samples <- expand_grid(Num = c(3L,6L,9L), Seed = names(cdata))
pb <- progress_bar$new(total = nrow(cdata_samples))
set.seed(26)
cdata_samples <- cdata_samples %>%
  mutate(Sample = map2(Num, Seed, ~ {pb$tick(); choose_complexes(.x, .y)}))

################################################################################################
## Results for HDBSCAN                                                                        ##
################################################################################################

MIN_PTS = c(3,5,10)

print("HDBSCAN, Euclidean distance")
pb <- progress_bar$new(total = nrow(cdata_samples))
results <- cdata_samples %>%
  mutate(Eucl = map(Sample, ~ {
    pb$tick();
    ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
    dist_obj <- dist(ds);
    d <- ncol(ds);
    run_hdbscan(dist_obj, d, MIN_PTS)
  }))

print("HDBSCAN, Pearson correlation")
pb <- progress_bar$new(total = nrow(results))
results <- results %>%
  mutate(Pear = map(Sample, ~ {
    pb$tick();
    ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
    pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
    dist_obj <- dist(ds, pear_dissim);
    d <- ncol(ds);
    run_hdbscan(dist_obj, d, MIN_PTS)
  }))

print("HDBSCAN, Parametric Euclidean Distance")
pb <- progress_bar$new(total = nrow(results))
results <- results %>%
  mutate(Eucl_params = map(Sample, ~ {
    pb$tick();
    ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
    ds <- make_unit_var(ds, c(Param_b, Param_c, Param_e));
    dist_obj <- dist(ds);
    d <- ncol(ds);
    run_hdbscan(dist_obj, d, MIN_PTS)
  }))

write_rds(results, here("thesis-scripts", "results", "hdbscan.rds"))


# pb <- progress_bar$new(total = 15) #
# results <- results %>%
#   mutate(Pear = map(Sample, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
#     dist_obj <- dist(ds, pear_dissim);
#     d <- ncol(ds);
#     run_hdbscan(dist_obj, d, MIN_PTS)
#   }))
# 
# pb <- progress_bar$new(total = nrow(cdata_samples))
# pb1 <- progress_bar$new(total = nrow(cdata_samples))
# pb2 <- progress_bar$new(total = nrow(cdata_samples))
# print("HDBSCAN clusterings, minPts = 3")
# results <- cdata_samples %>%
#   mutate(HDBSCAN3 = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     dist_obj <- dist(ds);
#     d <- 9;
#     min_pts <- 3;
#     run_hdbscan(dist_obj, d, min_pts)
#   })) %>%
#   mutate(HDBSCAN3_cor = map2(Sample, Num, ~ {
#     pb1$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
#     dist_obj <- dist(ds, pear_dissim);
#     d <- 9;
#     min_pts <- 3;
#     run_hdbscan(dist_obj, d, min_pts)
#   })) %>%
#   mutate(HDBSCAN3_p = map2(Sample, Num, ~ {
#     pb2$tick();
#     ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
#     ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e));
#     dist_obj <- dist(ds);
#     d <- 3;
#     min_pts <- 3;
#     run_hdbscan(dist_obj, d, min_pts)
#   }))
# 
# pb <- progress_bar$new(total = nrow(results))
# pb1 <- progress_bar$new(total = nrow(results))
# pb2 <- progress_bar$new(total = nrow(results))
# print("HDBSCAN clusterings, minPts = 5")
# results <- results %>%
#   mutate(HDBSCAN5 = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     dist_obj <- dist(ds);
#     d <- 9;
#     min_pts <- 5;
#     run_hdbscan(dist_obj, d, min_pts)
#   })) %>%
#   mutate(HDBSCAN5_cor = map2(Sample, Num, ~ {
#     pb1$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
#     dist_obj <- dist(ds, pear_dissim);
#     d <- 9;
#     min_pts <- 5;
#     run_hdbscan(dist_obj, d, min_pts)
#   })) %>%
#   mutate(HDBSCAN5_p = map2(Sample, Num, ~ {
#     pb2$tick();
#     ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
#     ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e));
#     dist_obj <- dist(ds);
#     d <- 3;
#     min_pts <- 5;
#     run_hdbscan(dist_obj, d, min_pts)
#   }))
# 
# pb <- progress_bar$new(total = nrow(results))
# pb1 <- progress_bar$new(total = nrow(results))
# pb2 <- progress_bar$new(total = nrow(results))
# print("HDBSCAN clusterings, minPts = 10")
# results <- results %>%
#   mutate(HDBSCAN10 = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     dist_obj <- dist(ds);
#     d <- 9;
#     min_pts <- 10;
#     run_hdbscan(dist_obj, d, min_pts)
#   })) %>%
#   mutate(HDBSCAN10_cor = map2(Sample, Num, ~ {
#     pb1$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
#     dist_obj <- dist(ds, pear_dissim);
#     d <- 9;
#     min_pts <- 10;
#     run_hdbscan(dist_obj, d, min_pts)
#   })) %>%
#   mutate(HDBSCAN10_p = map2(Sample, Num, ~ {
#     pb2$tick();
#     ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
#     ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e));
#     dist_obj <- dist(ds);
#     d <- 3;
#     min_pts <- 10;
#     run_hdbscan(dist_obj, d, min_pts)
#   }))
# 
# results <- results %>%
#   mutate(HDBSCAN_FOSC = pmap(list(HDBSCAN3, HDBSCAN5, HDBSCAN10), ~ {
#     res1 <- ..1$FOSC;
#     res2 <- ..2$FOSC;
#     res3 <- ..3$FOSC;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       best_cl = 1
#     }
#     res[[best_cl]]
#   })) %>%
#   mutate(HDBSCAN_cor_FOSC = pmap(list(HDBSCAN3_cor, HDBSCAN5_cor, HDBSCAN10_cor), ~ {
#     res1 <- ..1$FOSC;
#     res2 <- ..2$FOSC;
#     res3 <- ..3$FOSC;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       best_cl = 1
#     }
#     res[[best_cl]]
#   })) %>% 
#   mutate(HDBSCAN_p_FOSC = pmap(list(HDBSCAN3_p, HDBSCAN5_p, HDBSCAN10_p), ~ {
#     res1 <- ..1$FOSC;
#     res2 <- ..2$FOSC;
#     res3 <- ..3$FOSC;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       best_cl = 1
#     }
#     res[[best_cl]]
#   })) %>%
#   mutate(HDBSCAN_DBCV = pmap(list(HDBSCAN3, HDBSCAN5, HDBSCAN10, HDBSCAN_FOSC), ~ {
#     res1 <- ..1$HORIZ;
#     res2 <- ..2$HORIZ;
#     res3 <- ..3$HORIZ;
#     res4 <- ..4;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       res4
#     } else {
#       res[[best_cl]]
#     }
#   })) %>%
#   mutate(HDBSCAN_cor_DBCV = pmap(list(HDBSCAN3_cor, HDBSCAN5_cor, HDBSCAN10_cor, HDBSCAN_cor_FOSC), ~ {
#     res1 <- ..1$HORIZ;
#     res2 <- ..2$HORIZ;
#     res3 <- ..3$HORIZ;
#     res4 <- ..4;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       res4
#     } else {
#       res[[best_cl]]
#     }
#   })) %>%
#   mutate(HDBSCAN_p_DBCV = pmap(list(HDBSCAN3_p, HDBSCAN5_p, HDBSCAN10_p, HDBSCAN_p_FOSC), ~ {
#     res1 <- ..1$HORIZ;
#     res2 <- ..2$HORIZ;
#     res3 <- ..3$HORIZ;
#     res4 <- ..4;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       res4
#     } else {
#       res[[best_cl]]
#     }
#   })) %>%
#   select(Num, Seed, Sample, starts_with("HDBSCAN_"))
# 
# write_rds(results, here("thesis-scripts", "results", "hdbscan.rds"))
# 
# ################################################################################################
# ## Results for Single Linkage                                                                 ##
# ################################################################################################
# 
# set.seed(27)
# print("Single linkage clusterings")
# pb <- progress_bar$new(total = nrow(cdata_samples))
# pb1 <- progress_bar$new(total = nrow(cdata_samples))
# pb2 <- progress_bar$new(total = nrow(cdata_samples))
# print("Single linkage, minPts = 3")
# results1 <- cdata_samples %>%
#   mutate(Single3 = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     dist_obj <- dist(ds);
#     d <- 9;
#     min_pts <- 3;
#     run_single_linkage(dist_obj, d, min_pts)
#   })) %>%
#   mutate(Single3_cor = map2(Sample, Num, ~ {
#     pb1$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
#     dist_obj <- dist(ds, pear_dissim);
#     d <- 9;
#     min_pts <- 3;
#     run_single_linkage(dist_obj, d, min_pts)
#   })) %>%
#   mutate(Single3_p = map2(Sample, Num, ~ {
#     pb2$tick();
#     ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
#     ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e));
#     dist_obj <- dist(ds);
#     d <- 3;
#     min_pts <- 3;
#     run_single_linkage(dist_obj, d, min_pts)
#   }))
# 
# pb <- progress_bar$new(total = nrow(results1))
# pb1 <- progress_bar$new(total = nrow(results1))
# pb2 <- progress_bar$new(total = nrow(results1))
# print("Single linkage, minPts = 5")
# results1 <- results1 %>%
#   mutate(Single5 = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     dist_obj <- dist(ds);
#     d <- 9;
#     min_pts <- 5;
#     run_single_linkage(dist_obj, d, min_pts)
#   })) %>%
#   mutate(Single5_cor = map2(Sample, Num, ~ {
#     pb1$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
#     dist_obj <- dist(ds, pear_dissim);
#     d <- 9;
#     min_pts <- 5;
#     run_single_linkage(dist_obj, d, min_pts)
#   })) %>%
#   mutate(Single5_p = map2(Sample, Num, ~ {
#     pb2$tick();
#     ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
#     ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e));
#     dist_obj <- dist(ds);
#     d <- 3;
#     min_pts <- 5;
#     run_single_linkage(dist_obj, d, min_pts)
#   }))
# 
# pb <- progress_bar$new(total = nrow(results1))
# pb1 <- progress_bar$new(total = nrow(results1))
# pb2 <- progress_bar$new(total = nrow(results1))
# print("Single linkage, minPts = 10")
# results1 <- results1 %>%
#   mutate(Single10 = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     dist_obj <- dist(ds);
#     d <- 9;
#     min_pts <- 10;
#     run_single_linkage(dist_obj, d, min_pts)
#   })) %>%
#   mutate(Single10_cor = map2(Sample, Num, ~ {
#     pb1$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     pear_dissim <- function(x1, x2) (1 - cor(x1, x2)) / 2;
#     dist_obj <- dist(ds, pear_dissim);
#     d <- 9;
#     min_pts <- 10;
#     run_single_linkage(dist_obj, d, min_pts)
#   })) %>%
#   mutate(Single10_p = map2(Sample, Num, ~ {
#     pb2$tick();
#     ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
#     ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e));
#     dist_obj <- dist(ds);
#     d <- 3;
#     min_pts <- 10;
#     run_single_linkage(dist_obj, d, min_pts)
#   }))
# 
# results1 <- results1 %>%
#   mutate(Single_FOSC = pmap(list(Single3, Single5, Single10), ~ {
#     res1 <- ..1$FOSC;
#     res2 <- ..2$FOSC;
#     res3 <- ..3$FOSC;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       best_cl = 1
#     }
#     res[[best_cl]]
#   })) %>%
#   mutate(Single_cor_FOSC = pmap(list(Single3_cor, Single5_cor, Single10_cor), ~ {
#     res1 <- ..1$FOSC;
#     res2 <- ..2$FOSC;
#     res3 <- ..3$FOSC;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       best_cl = 1
#     }
#     res[[best_cl]]
#   })) %>% 
#   mutate(Single_p_FOSC = pmap(list(Single3_p, Single5_p, Single10_p), ~ {
#     res1 <- ..1$FOSC;
#     res2 <- ..2$FOSC;
#     res3 <- ..3$FOSC;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       best_cl = 1
#     }
#     res[[best_cl]]
#   })) %>%
#   mutate(Single_DBCV = pmap(list(Single3, Single5, Single10, Single_FOSC), ~ {
#     res1 <- ..1$HORIZ;
#     res2 <- ..2$HORIZ;
#     res3 <- ..3$HORIZ;
#     res4 <- ..4;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       res4
#     } else {
#       res[[best_cl]]
#     }
#   })) %>%
#   mutate(Single_cor_DBCV = pmap(list(Single3_cor, Single5_cor, Single10_cor, Single_cor_FOSC), ~ {
#     res1 <- ..1$HORIZ;
#     res2 <- ..2$HORIZ;
#     res3 <- ..3$HORIZ;
#     res4 <- ..4;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       res4
#     } else {
#       res[[best_cl]]
#     }
#   })) %>%
#   mutate(Single_p_DBCV = pmap(list(Single3_p, Single5_p, Single10_p, Single_p_FOSC), ~ {
#     res1 <- ..1$HORIZ;
#     res2 <- ..2$HORIZ;
#     res3 <- ..3$HORIZ;
#     res4 <- ..4;
#     res <- list(res1, res2, res3);
#     dbcv <- map_dbl(res, ~ .x$dbcv);
#     best_cl <- which.max(dbcv);
#     if (length(best_cl) == 0) {
#       res4
#     } else {
#       res[[best_cl]]
#     }
#   })) %>% select(Num, Seed, Sample, starts_with("Single_"))
# 
# write_rds(results1, here("thesis-scripts", "results", "single_linkage.rds"))
# 
# 
# i = 535
# ds <- cdata_samples$Sample[[i]]
# ds <- (cdata_samples %>% filter(Num == 3 & Seed == "3040"))$Sample[[1]]
# ds
# ds <- ds %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37)
# dist_obj <- dist(ds)
# d <- ncol(ds)
# min_pts <- c(3, 5, 10)
# res <- run_hdbscan(dist_obj, d, min_pts)
# res$single_cut
# res$hdbscan_cut
# res$cut
# res$fosc


# i = 4
# set.seed(27)
# ds <- cdata_samples %>% slice_sample(n = 10) %>% pull(Sample) %>% pluck(i)
# ds
# ds <- ds %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37)
# # ds <- ds %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein)
# # ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e))
# dist_obj <- dist(ds)
# d <- 9
# # d <- 3
# min_pts <- 3
# res <- run_hdbscan(dist_obj, d, min_pts)
# res1 <- run_single_linkage(dist_obj, d, min_pts)
# 
# i = 4
# ds <- results$Sample[[i]]
# ds <- ds %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37)
# dist_obj <- dist(ds)
# d <- 9
# min_pts <- 3
# res <- run_single_linkage(dist_obj, d, min_pts)
