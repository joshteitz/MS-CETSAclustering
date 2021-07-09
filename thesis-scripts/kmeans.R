library(conflicted)
library(tidyverse)
library(here)
library(progress)

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
## Results for K-means                                                                        ##
################################################################################################

# For each sample, run k-means.
print("Run K-means for each dataset")
set.seed(27)
pb <- progress_bar$new(total = nrow(cdata_samples))
pb1 <- progress_bar$new(total = nrow(cdata_samples))
results <- cdata_samples %>%
  mutate(Kmeans = map2(Sample, Num, ~ {
    pb$tick();
    ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein);
    run_kmeans(ds, .y * 2)
  })) %>%
  mutate(Kmeans_p = map2(Sample, Num, ~ {
    pb1$tick();
    ds <- .x %>% inner_join(mparams, by = "Protein") %>% select(-Complex, -Protein);
    ds <- make_unit_var(ds, numeric_cols = c(Param_b, Param_c, Param_e))
    run_kmeans(ds, .y * 2)
  })) %>%
  mutate(ARI = map2_dbl(Sample, Kmeans, ~ {
    gt <- .x$Complex;
    lbls <- .y$cl;
    adjustedRandIndex(gt, lbls)
  })) %>%
  mutate(ARI_p = map2_dbl(Sample, Kmeans_p, ~ {
    gt <- .x$Complex;
    lbls <- .y$cl;
    adjustedRandIndex(gt, lbls)
  }))

# Save the results
write_rds(results, here("thesis-scripts", "results", "kmeans.rds"))

# # For each clustering, calculate the adjusted RAND index and other external evaluations
# pb <- progress_bar$new(total = nrow(results))
# results <- results %>%
#   mutate(ARI = map2_dbl(Sample, Kmeans, ~ get_ARI(.x$Complex, .y$cl))) %>%
#   mutate(Extern_eval = map2(Sample, Kmeans, ~ {pb$tick(); external_eval(.x$Complex, .y$cl)}))
# 
# # For each cluster, calculate its average silhouette width and maximum F-measure
# cluster_scores <- results %>%
#   mutate(Cluster_Avg_sil = map(Kmeans, ~ {
#     .x$sils %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Cluster_Fmeas = map(Extern_eval, ~ .x$maxfmeas)) %>%
#   select(-Sample, - Kmeans, -ARI, -Extern_eval) %>%
#   mutate(Cluster_Avg_sil_Fmeas = map2(Cluster_Avg_sil, Cluster_Fmeas, ~ {
#     .x %>% 
#       inner_join(.y, by = "Cluster") %>%
#       distinct(Cluster, .keep_all = T)
#   })) %>%
#   select(-Cluster_Avg_sil, -Cluster_Fmeas) %>%
#   unnest(Cluster_Avg_sil_Fmeas)
# 
# write_rds(results, paste0(here("thesis-scripts", "data", "kmeans"), "/results.rds"))
# write_rds(cluster_scores, paste0(here("thesis-scripts", "data", "kmeans"), "/cluster_scores.rds"))
# 
# results <- read_rds(paste0(here("thesis-scripts", "data", "kmeans"), "/results.rds"))
# cluster_scores <- read_rds(paste0(here("thesis-scripts", "data", "kmeans"), "/cluster_scores.rds"))
# 
# ################################################################################################
# ## Results for K-means on melting parameter data                                              ##
# ################################################################################################
# 
# # For each sample, run k-means.
# pb <- progress_bar$new(total = nrow(cdata_samples))
# results1 <- cdata_samples %>%
#   mutate(Kmeans = map2(Sample, Num, ~ {
#     pb$tick();
#     ds <- .x %>% inner_join(mdata, by = "Protein") %>% select(-Complex, -Protein);
#     run_kmeans(ds, .y * 2)
#   }))
# 
# # For each clustering, calculate the adjusted RAND index and other external evaluations
# pb <- progress_bar$new(total = nrow(results1))
# results1 <- results1 %>%
#   mutate(ARI = map2_dbl(Sample, Kmeans, ~ get_ARI(.x$Complex, .y$cl))) %>%
#   mutate(Extern_eval = map2(Sample, Kmeans, ~ {pb$tick(); external_eval(.x$Complex, .y$cl)}))
# 
# # For each cluster, calculate its average silhouette width and maximum F-measure
# cluster_scores1 <- results1 %>%
#   mutate(Cluster_Avg_sil = map(Kmeans, ~ {
#     .x$sils %>%
#       group_by(Cluster) %>%
#       summarize(Avg_sil = mean(Sil), .groups = "drop") %>%
#       arrange(Cluster)
#   })) %>%
#   mutate(Cluster_Fmeas = map(Extern_eval, ~ .x$maxfmeas)) %>%
#   select(-Sample, - Kmeans, -ARI, -Extern_eval) %>%
#   mutate(Cluster_Avg_sil_Fmeas = map2(Cluster_Avg_sil, Cluster_Fmeas, ~ {
#     .x %>% 
#       inner_join(.y, by = "Cluster") %>%
#       distinct(Cluster, .keep_all = T)
#   })) %>%
#   select(-Cluster_Avg_sil, -Cluster_Fmeas) %>%
#   unnest(Cluster_Avg_sil_Fmeas)
# 
# write_rds(results1, paste0(here("thesis-scripts", "data", "kmeans"), "/results1.rds"))
# write_rds(cluster_scores1, paste0(here("thesis-scripts", "data", "kmeans"), "/cluster_scores1.rds"))
# 
# results1 <- read_rds(paste0(here("thesis-scripts", "data", "kmeans"), "/results1.rds"))
# cluster_scores1 <- read_rds(paste0(here("thesis-scripts", "data", "kmeans"), "/cluster_scores1.rds"))
# 
# ################################################################################################
# ## Visualization of Results                                                                   ##
# ################################################################################################
# 
# CLUSTER_SCORES <- cluster_scores
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
# # Correlation of Average silhouette width vs. F-measure
# NUM = 3
# cor(
#   CLUSTER_SCORES %>% filter(Num == NUM) %>% pull(Avg_sil),
#   CLUSTER_SCORES %>% filter(Num == NUM) %>% pull(Fmeas)
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

#################################################################################

# # Data for box plots
# plot_data <- results1 %>%
#   mutate(MaxF = map_dbl(Extern_eval, ~ .x$maxfmeas$Fmeas %>% max)) %>%
#   select(Seed, Num, ARI, MaxF) %>%
#   mutate(Num = Num %>% as.character)
# 
# # ARI box plot
# ggplot(plot_data, aes(x = Num, y = ARI)) +
#   geom_boxplot(aes(fill = Num), outlier.size = 3, lwd = 1.5, fatten = 1) +
#   scale_y_continuous(breaks = seq(-1, 1, length.out = 11)) +
#   # expand_limits(y = c(-.9, .9)) +
#   xlab("Datasets") +
#   scale_x_discrete(labels = c("3" = "3C", "6" = "6C", "9" = "9C")) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 28),
#         axis.title = element_text(size = 28),
#         legend.position = "none") +
#   scale_fill_brewer(palette = "Blues")
# 
# # Maximum F-measure box plot
# ggplot(plot_data, aes(x = Num, y = MaxF)) +
#   geom_boxplot(aes(fill = Num), outlier.size = 3, lwd = 1.5, fatten = 1) +
#   scale_y_continuous(breaks = seq(-1, 1, length.out = 11)) +
#   # expand_limits(y = c(-.9, .9)) +
#   xlab("Datasets") +
#   ylab("Maximum F-measure") +
#   scale_x_discrete(labels = c("3" = "3C", "6" = "6C", "9" = "9C")) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 28),
#         axis.title = element_text(size = 28),
#         legend.position = "none") +
#   scale_fill_brewer(palette = "Blues")
# 
# # scatter plot of Average silhouette width vs. F-measure. Each point represents a single cluster.
# NUM = 9
# ggplot(cluster_scores %>% filter(Num == NUM), aes(x = Avg_sil, y = Fmeas)) +
#   geom_point(size = 2) +
#   scale_y_continuous(breaks = seq(-1, 1, length.out = 21)) +
#   scale_x_continuous(breaks = seq(-1, 1, length.out = 21)) +
#   expand_limits(y = c(0,1), x = c(0,1)) +
#   theme_bw() +
#   theme(axis.text = element_text(size = 28),
#         axis.title = element_text(size = 28))

# NUM = 9
# cor(cluster_scores %>% filter(Num == NUM) %>% pull(Avg_sil), cluster_scores %>% filter(Num == NUM) %>% pull(Fmeas))
# cor(cluster_scores %>% pull(Avg_sil), cluster_scores %>% pull(Fmeas))
# 
# # Label a cluster as 1 if Fmeas >= 2/3 and 0 if Fmeas < 2/3.
# cluster_scores <- cluster_scores %>%
#   mutate(Label = map_int(Fmeas, ~ ifelse(.x >= 2/3, 1L, 0L)))
# 
# par(pty = "s")
# roc_plot <- roc((cluster_scores)$Label, 
#     (cluster_scores)$Avg_sil, 
#     plot = T, 
#     legacy.axes = T,
#     xlab = "False Positive Rate",
#     ylab = "True Positive Rate",
#     col = "#477eb8", lwd = 8,
#     print.auc = T,
#     cex.lab = 2, cex.axis = 2)
# 
# roc_info <- tibble(
#   TPR = roc_plot$sensitivities,
#   FPR = (1 - roc_plot$specificities) * 100
# )