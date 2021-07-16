library(conflicted)
library(tidyverse)
library(here)
library(progress)
library(pROC)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

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
## Results for GMMs                                                                           ##
################################################################################################
set.seed(27)
print("Fit GMMs to each dataset")
pb <- progress_bar$new(total = nrow(cdata_samples))
pb1 <- progress_bar$new(total = nrow(cdata_samples))
results <- cdata_samples %>%
  mutate(GMM = map2(Sample, Num, ~ {
    pb$tick();
    ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
    run_gmms(ds, .y * 2)
  })) %>%
  mutate(GMM_p = map2(Sample, Num, ~ {
    pb1$tick();
    ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
    ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e));
    run_gmms(ds, .y * 2)
  }))

print("Compute silhouettes w.r.t Euclidean distance")
pb <- progress_bar$new(total = nrow(results))
results <- results %>%
  mutate(Sils = map2(Sample, GMM, ~ {
    pb$tick();
    lbls <- .y$cl;
    ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
    dist_obj <- dist(ds);
    get_sils(lbls, dist_obj)
  }))

print("Compute silhouettes w.r.t Parametric Euclidean distance")
pb <- progress_bar$new(total = nrow(results))
results <- results %>%
  mutate(Sils_p = map2(Sample, GMM_p, ~ {
    pb$tick();
    lbls <- .y$cl;
    ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
    ds <- make_unit_var(ds, c(Param_b, Param_c, Param_e));
    dist_obj <- dist(ds);
    get_sils(lbls, dist_obj)
  }))

print("Compute silhouettes w.r.t Mahalanobis distance")
pb <- progress_bar$new(total = nrow(results))
results <- results %>%
  mutate(Sils_mahal = map2(Sample, GMM, ~{
    pb$tick();
    lbls <- .y$cl;
    ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
    sigmas <- .y$sigmas;
    get_sils_mahal(lbls, ds, sigmas)
  }))

print("Compute silhouettes w.r.t. Parametric Mahalanobis distance")
pb <- progress_bar$new(total = nrow(results))
results <- results %>%
  mutate(Sils_mahal_p = map2(Sample, GMM_p, ~{
    pb$tick();
    lbls <- .y$cl;
    ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
    ds <- make_unit_var(ds, c(Param_b, Param_c, Param_e));
    sigmas <- .y$sigmas;
    get_sils_mahal(lbls, ds, sigmas)
  }))

print("Compute ARI for each dataset w.r.t. the ground truth")
results <- results %>%
  mutate(ARI = map2_dbl(Sample, GMM, ~ {
    gt <- .x$Complex;
    cl <- .y$cl;
    adjustedRandIndex(gt, cl)
  })) %>%
  mutate(ARI_p = map2_dbl(Sample, GMM_p, ~ {
    gt <- .x$Complex;
    cl <- .y$cl;
    adjustedRandIndex(gt, cl)
  }))

write_rds(results, here("thesis-scripts", "results", "gmm.rds"))

# # For each sample, run GMMs.
# pb <- progress_bar$new(total = nrow(cdata_samples))
# results <- cdata_samples %>% slice_sample(n = 10) %>%
#   mutate(GMM = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein, -T37);
#     run_gmms(ds, .y * 2)
#   }))
# 
# # Add a logical column called `ICL`. The column is TRUE iff ICL selects a different model than BIC.
# results <- results %>%
#   mutate(ICL = map_lgl(GMM, ~ !is.null(.x$ICL)))
# 
# # For each clustering, calculate its silhouette values, first using Euclidean distance as the
# # dissimilarity measure, and then using Mahalanobis as the dissimilarity measure.
# pb1 <- progress_bar$new(total = nrow(results))
# pb2 <- progress_bar$new(total = nrow(results))
# results <- results %>%
#   mutate(Sils = map2(GMM, Sample, ~ {
#     pb1$tick();
#     lbls <- .x$BIC$cl;
#     ds <- .y %>% inner_join(mdata, by = "Protein") %>% select(-Protein, -Complex, -T37);
#     get_sils(lbls, dist(ds))
#   })) %>%
#   mutate(Sils_mahal = map2(GMM, Sample, ~ {
#     pb2$tick();
#     lbls <- .x$BIC$cl;
#     ds <- .y %>% inner_join(mdata, by = "Protein") %>% select(-Protein, -Complex, -T37);
#     sigmas <- .x$BIC$sigmas;
#     get_sils_mahal(lbls, ds, sigmas)
#   }))
# 
# # For each clustering, calculate the adjusted RAND index and other external evaluations
# pb <- progress_bar$new(total = nrow(results))
# results <- results %>%
#   mutate(ARI = map2_dbl(Sample, GMM, ~ get_ARI(.x$Complex, .y$BIC$cl))) %>%
#   mutate(Extern_eval = map2(Sample, GMM, ~ {pb$tick(); external_eval(.x$Complex, .y$BIC$cl)}))
# 
# # For each cluster, calculate its average silhouette width and maximum F-measure
# cluster_scores <- results %>%
#   mutate(Avg_sil = map(Sils, ~ {
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Avg_sil_mahal = map(Sils_mahal, ~ {
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil_mahal = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Fmeas = map(Extern_eval, ~ .x$maxfmeas)) %>%
#   select(-Sample, - GMM, -ICL, -Sils, -Sils_mahal, -ARI, -Extern_eval) %>%
#   mutate(Scores = pmap(list(Avg_sil, Avg_sil_mahal, Fmeas), ~ {
#     bind_cols(
#       ..1,
#       ..2 %>% select(Avg_sil_mahal),
#       ..3 %>% distinct(Cluster, .keep_all = T) %>% select(Fmeas)
#     )
#   })) %>%
#   select(-Avg_sil, -Avg_sil_mahal, -Fmeas) %>%
#   unnest(Scores)
# 
# write_rds(results, paste0(here("thesis-scripts", "data", "gmm"), "/results.rds"))
# write_rds(cluster_scores, paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores.rds"))
# 
# results <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/results.rds"))
# cluster_scores <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores.rds"))
# 
# ################################################################################################
# ## Results of GMMs on melting curve data - ICL selection                                      ##
# ################################################################################################
# 
# results1 <- results
# 
# # For each clustering, calculate its silhouette values, first using Euclidean distance as the
# # dissimilarity measure, and then using Mahalanobis as the dissimilarity measure. 
# pb1 <- progress_bar$new(total = nrow(results1))
# pb2 <- progress_bar$new(total = nrow(results1))
# results1 <- results1 %>%
#   mutate(Sils = pmap(list(GMM, Sample, ICL, Sils), ~ {
#     pb1$tick();
#     if (..3) {
#       lbls <- ..1$ICL$cl;
#       ds <- ..2 %>% inner_join(mdata, by = "Protein") %>% select(-Protein, -Complex, -T37);
#       get_sils(lbls, dist(ds))
#     } else {
#       ..4
#     }
#   })) %>%
#   mutate(Sils_mahal = pmap(list(GMM, Sample, ICL, Sils_mahal), ~ {
#     pb2$tick();
#     if (..3) {
#       lbls <- ..1$ICL$cl;
#       ds <- ..2 %>% inner_join(mdata, by = "Protein") %>% select(-Protein, -Complex, -T37);
#       sigmas <- ..1$ICL$sigmas;
#       get_sils_mahal(lbls, ds, sigmas)
#     } else {
#       ..4
#     }
#   }))
# 
# # For each clustering, calculate the adjusted RAND index and other external evaluations
# pb <- progress_bar$new(total = nrow(results1))
# results1 <- results1 %>%
#   mutate(ARI = pmap_dbl(list(Sample, GMM, ICL, ARI), ~ {
#     if (..3) {
#       get_ARI(..1$Complex, ..2$ICL$cl)
#     } else {
#       ..4
#     }
#   })) %>%
#   mutate(Extern_eval = pmap(list(Sample, GMM, ICL, Extern_eval), ~ {
#     pb$tick();
#     if (..3) {
#       external_eval(..1$Complex, ..2$ICL$cl)
#     } else {
#       ..4
#     }
#   }))
# 
# # For each cluster, calculate its average silhouette width and maximum F-measure
# pb1 <- progress_bar$new(total = nrow(results1))
# pb2 <- progress_bar$new(total = nrow(results1))
# cluster_scores1 <- results1 %>%
#   mutate(Avg_sil = map(Sils, ~ {
#     pb1$tick();
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Avg_sil_mahal = map(Sils_mahal, ~ {
#     pb2$tick();
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil_mahal = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Fmeas = map(Extern_eval, ~ .x$maxfmeas)) %>%
#   select(-Sample, - GMM, -ICL, -Sils, -Sils_mahal, -ARI, -Extern_eval) %>%
#   mutate(Scores = pmap(list(Avg_sil, Avg_sil_mahal, Fmeas), ~ {
#     bind_cols(
#       ..1,
#       ..2 %>% select(Avg_sil_mahal),
#       ..3 %>% distinct(Cluster, .keep_all = T) %>% select(Fmeas)
#     )
#   })) %>%
#   select(-Avg_sil, -Avg_sil_mahal, -Fmeas) %>%
#   unnest(Scores)
# 
# write_rds(results1, paste0(here("thesis-scripts", "data", "gmm"), "/results1.rds"))
# write_rds(cluster_scores1, paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores1.rds"))
# 
# results1 <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/results1.rds"))
# cluster_scores1 <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores1.rds"))
# 
# ################################################################################################
# ## Results of GMMs on melting parameter data - BIC selection                                  ##
# ################################################################################################
# 
# # For each sample, run GMMs.
# pb <- progress_bar$new(total = nrow(cdata_samples))
# results2 <- cdata_samples %>% slice_sample(n = 10) %>%
#   mutate(GMM = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
#     run_gmms(ds, .y * 2)
#   }))
# 
# # Add a logical column called `ICL`. The column is TRUE iff ICL selects a different model than BIC.
# results2 <- results2 %>%
#   mutate(ICL = map_lgl(GMM, ~ !is.null(.x$ICL)))
# 
# # For each clustering, calculate its silhouette values, first using Euclidean distance as the
# # dissimilarity measure, and then using Mahalanobis as the dissimilarity measure.
# pb1 <- progress_bar$new(total = nrow(results2))
# pb2 <- progress_bar$new(total = nrow(results2))
# results2 <-
#   results2 %>%
#   mutate(Sils = map2(GMM, Sample, ~ {
#     pb1$tick();
#     lbls <- .x$BIC$cl;
#     ds <- .y %>% inner_join(mparams, by = "Protein") %>% select(-Protein, -Complex);
#     get_sils(lbls, dist(ds))
#   })) %>%
#   mutate(Sils_mahal = map2(GMM, Sample, ~ {
#     pb2$tick();
#     lbls <- .x$BIC$cl;
#     ds <- .y %>% inner_join(mparams, by = "Protein") %>% select(-Protein, -Complex);
#     sigmas <- .x$BIC$sigmas;
#     get_sils_mahal(lbls, ds, sigmas)
#   }))
# 
# # For each clustering, calculate the adjusted RAND index and other external evaluations
# pb <- progress_bar$new(total = nrow(results2))
# results2 <- results2 %>%
#   mutate(ARI = map2_dbl(Sample, GMM, ~ get_ARI(.x$Complex, .y$BIC$cl))) %>%
#   mutate(Extern_eval = map2(Sample, GMM, ~ {pb$tick(); external_eval(.x$Complex, .y$BIC$cl)}))
# 
# # For each cluster, calculate its average silhouette width and maximum F-measure
# pb1 <- progress_bar$new(total = nrow(results2))
# pb2 <- progress_bar$new(total = nrow(results2))
# cluster_scores2 <- results2 %>%
#   mutate(Avg_sil = map(Sils, ~ {
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Avg_sil_mahal = map(Sils_mahal, ~ {
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil_mahal = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Fmeas = map(Extern_eval, ~ .x$maxfmeas)) %>%
#   select(-Sample, - GMM, -ICL, -Sils, -Sils_mahal, -ARI, -Extern_eval) %>%
#   mutate(Scores = pmap(list(Avg_sil, Avg_sil_mahal, Fmeas), ~ {
#     bind_cols(
#       ..1,
#       ..2 %>% select(Avg_sil_mahal),
#       ..3 %>% distinct(Cluster, .keep_all = T) %>% select(Fmeas)
#     )
#   })) %>%
#   select(-Avg_sil, -Avg_sil_mahal, -Fmeas) %>%
#   unnest(Scores)
# 
# write_rds(results2, paste0(here("thesis-scripts", "data", "gmm"), "/results2.rds"))
# write_rds(cluster_scores2, paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores2.rds"))
# 
# results2 <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/results2.rds"))
# cluster_scores2 <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores2.rds"))
# 
# ################################################################################################
# ## Results of GMMs on melting parameter data - ICL selection                                  ##
# ################################################################################################
# 
# results3 <- results2
# 
# # For each clustering, calculate its silhouette values, first using Euclidean distance as the
# # dissimilarity measure, and then using Mahalanobis as the dissimilarity measure. 
# pb1 <- progress_bar$new(total = nrow(results3))
# pb2 <- progress_bar$new(total = nrow(results3))
# results3 <- results3 %>%
#   mutate(Sils = pmap(list(GMM, Sample, ICL, Sils), ~ {
#     pb1$tick();
#     if (..3) {
#       lbls <- ..1$ICL$cl;
#       ds <- ..2 %>% inner_join(mparams, by = "Protein") %>% select(-Protein, -Complex);
#       get_sils(lbls, dist(ds))
#     } else {
#       ..4
#     }
#   })) %>%
#   mutate(Sils_mahal = pmap(list(GMM, Sample, ICL, Sils_mahal), ~ {
#     pb2$tick();
#     if (..3) {
#       lbls <- ..1$ICL$cl;
#       ds <- ..2 %>% inner_join(mparams, by = "Protein") %>% select(-Protein, -Complex);
#       sigmas <- ..1$ICL$sigmas;
#       get_sils_mahal(lbls, ds, sigmas)
#     } else {
#       ..4
#     }
#   }))
# 
# # For each clustering, calculate the adjusted RAND index and other external evaluations
# pb <- progress_bar$new(total = nrow(results3))
# results3 <- results3 %>%
#   mutate(ARI = pmap_dbl(list(Sample, GMM, ICL, ARI), ~ {
#     if (..3) {
#       get_ARI(..1$Complex, ..2$ICL$cl)
#     } else {
#       ..4
#     }
#   })) %>%
#   mutate(Extern_eval = pmap(list(Sample, GMM, ICL, Extern_eval), ~ {
#     pb$tick();
#     if (..3) {
#       external_eval(..1$Complex, ..2$ICL$cl)
#     } else {
#       ..4
#     }
#   }))
# 
# # For each cluster, calculate its average silhouette width and maximum F-measure
# pb1 <- progress_bar$new(total = nrow(results3))
# pb2 <- progress_bar$new(total = nrow(results3))
# cluster_scores3 <- results3 %>%
#   mutate(Avg_sil = map(Sils, ~ {
#     pb1$tick();
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Avg_sil_mahal = map(Sils_mahal, ~ {
#     pb2$tick();
#     .x %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil_mahal = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Fmeas = map(Extern_eval, ~ .x$maxfmeas)) %>%
#   select(-Sample, - GMM, -ICL, -Sils, -Sils_mahal, -ARI, -Extern_eval) %>%
#   mutate(Scores = pmap(list(Avg_sil, Avg_sil_mahal, Fmeas), ~ {
#     bind_cols(
#       ..1,
#       ..2 %>% select(Avg_sil_mahal),
#       ..3 %>% distinct(Cluster, .keep_all = T) %>% select(Fmeas)
#     )
#   })) %>%
#   select(-Avg_sil, -Avg_sil_mahal, -Fmeas) %>%
#   unnest(Scores)
# 
# write_rds(results3, paste0(here("thesis-scripts", "data", "gmm"), "/results3.rds"))
# write_rds(cluster_scores3, paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores3.rds"))
# 
# results3 <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/results3.rds"))
# cluster_scores3 <- read_rds(paste0(here("thesis-scripts", "data", "gmm"), "/cluster_scores3.rds"))
# 
# ################################################################################################
# ## Visualization of Results                                                                   ##
# ################################################################################################
# 
# CLUSTER_SCORES <- cluster_scores2
# 
# # scatter plot of Average silhouette width vs. F-measure. Each point represents a single cluster.
# num_labs <- c("3C", "6C", "9C")
# names(num_labs) <- c("3", "6", "9")
# ggplot(CLUSTER_SCORES, aes(x = Avg_sil, y = Fmeas)) +
#   geom_point(size = 2) +
#   geom_smooth(method = 'lm', formula = y ~ x, size = 2) +
#   scale_y_continuous(breaks = seq(-1, 1, length.out = 21)) +
#   scale_x_continuous(breaks = seq(-1, 1, length.out = 21)) +
#   expand_limits(y = c(0,1), x = c(0,1)) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         strip.text = element_text(size = 20)) +
#   facet_grid(cols = vars(Num), labeller = labeller(Num = num_labs))
# 
# # Scatter plot of Average silhouette width w.r.t. Mahalanobis distance vs. F-measure.
# # Each point represents a single cluster.
# num_labs <- c("3C", "6C", "9C")
# names(num_labs) <- c("3", "6", "9")
# ggplot(CLUSTER_SCORES, aes(x = Avg_sil_mahal, y = Fmeas)) +
#   geom_point(size = 2) +
#   geom_smooth(method = 'lm', formula = y ~ x, size = 2) +
#   scale_y_continuous(breaks = seq(-1, 1, length.out = 21)) +
#   scale_x_continuous(breaks = seq(-1, 1, length.out = 21)) +
#   expand_limits(y = c(0,1), x = c(0,1)) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         strip.text = element_text(size = 20)) +
#   facet_grid(cols = vars(Num), labeller = labeller(Num = num_labs))
# 
# # Correlation of Average silhouette width vs. F-measure
# NUM = 3
# cor(
#   CLUSTER_SCORES %>% filter(Num == NUM) %>% pull(Avg_sil),
#   CLUSTER_SCORES %>% filter(Num == NUM) %>% pull(Fmeas)
# )
# 
# # Correlation of Average silhouette width w.r.t. Mahalanobis distance
# # vs. F-measure
# NUM = 9
# cor(
#   CLUSTER_SCORES %>% filter(Num == NUM) %>% pull(Avg_sil_mahal),
#   cluster_scores %>% filter(Num == NUM) %>% pull(Fmeas)
# )
# 
# # Table of F-measures over different silhouette ranges.
# NUM = 3
# CLUSTER_SCORES %>%
#   filter(Num == NUM) %>%
#   group_by(Sil_rng = cut(Avg_sil, breaks = c(0,.25,.5,.7,1), right = F)) %>%
#   summarize(
#     n = length(Fmeas),
#     `F-meas < .5` = sum(Fmeas < .5) / length(Fmeas),
#     `F-meas >= .5` = sum(Fmeas >= .5) / length(Fmeas),
#     .groups = "drop"
#   )
# 
# # Table of F-measures over different silhouette ranges (Mahalanobis distance)
# NUM = 9
# CLUSTER_SCORES %>%
#   filter(Num == NUM) %>%
#   group_by(Sil_rng = cut(Avg_sil_mahal, breaks = c(0,.25,.5,.7,1), right = F)) %>%
#   summarize(
#     n = length(Fmeas),
#     `F-meas < .5` = sum(Fmeas < .5) / length(Fmeas),
#     `F-meas >= .5` = sum(Fmeas >= .5) / length(Fmeas),
#     .groups = "drop"
#   )
