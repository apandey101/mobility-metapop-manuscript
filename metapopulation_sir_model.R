###############################################################################
# Stochastic Metapopulation SIR Model for South Carolina Counties
#
# This script simulates disease spread across 46 SC counties using a
# stochastic SIR (Susceptible-Infected-Recovered) metapopulation framework.
# Between-patch transmission is modeled using two alternative approaches:
#
#   1. Radiation Model  — coupling derived from geographic distance and
#                         population via the classical radiation model
#   2. Commuting Model  — coupling derived from county-to-county commuting
#                         flow data
#
# For each approach, the script:
#   - Computes pairwise effective distances between counties
#   - Runs N stochastic simulations with multi-patch seeding
#   - Estimates the average first-case arrival time per county
#   - Correlates effective distance with arrival time (Spearman's rho)
#   - Produces labeled scatter plots of effective distance vs. onset
#
# Required input files (place in working directory or update paths):
#   - SC_county_sociodemographic.csv   : County populations & demographics
#   - sc_shape_files.RData             : SF geometry for SC counties
#   - commute_SC_clean.csv             : County-to-county commuting flows
#

###############################################################################

rm(list = ls())

# =============================================================================
# 1. LOAD PACKAGES
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(sf)
library(igraph)
library(future.apply)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

# County-level sociodemographic data (must contain 'County' and 'Total_population')
socio <- read.csv("SC_county_sociodemographic.csv")
socio <- socio %>% rename(county = County)

# Named population vector (used throughout)
N <- setNames(as.numeric(socio$Total_population), socio$county)

# SC county shapefiles (provides sf object: sc_county_sf)
load("sc_shape_files.RData")

# =============================================================================
# 3. GLOBAL SIMULATION PARAMETERS
# =============================================================================

n_patches     <- 46     # Number of SC counties
gamma         <- 0.25   # Recovery rate (1/gamma = 4-day infectious period)
n_simulations <- 1000   # Number of stochastic replicates

# =============================================================================
# 4. GEOGRAPHIC DISTANCE MATRIX
# =============================================================================
# Compute pairwise distances between county centroids (meters).
# Used by the radiation model to define spatial coupling.

counties_sf       <- st_transform(sc_county_sf, crs = 3857) # Project to Web Mercator
county_centroids  <- st_centroid(counties_sf)
dist_matrix       <- st_distance(county_centroids)

colnames(dist_matrix) <- rownames(dist_matrix) <- counties_sf$NAME

# Long-format distance table
distance_df <- as.data.frame(as.table(dist_matrix))
colnames(distance_df) <- c("county1", "county2", "distance_meters")
distance_df$distance_meters <- as.numeric(distance_df$distance_meters)

# =============================================================================
# 5. HELPER FUNCTIONS
# =============================================================================

#' Classical Radiation Model
#'
#' Computes pairwise coupling probabilities between patches using the
#' radiation model (Simini et al., 2012). For each pair (i, j), the
#' intervening opportunities s_ij are defined as the total population
#' of all counties closer to i than j is (excluding i and j).
#'
#' @param N          Named numeric vector of county populations
#' @param distance_df  Long-format data frame with county1, county2, distance_meters
#' @return Matrix of coupling probabilities (not yet normalised to transmission rates)

compute_radiation_classical <- function(N, distance_df) {
  counties  <- names(N)
  n         <- length(counties)
  beta_mat  <- matrix(0, n, n, dimnames = list(counties, counties))

  for (i in counties) {
    for (j in counties) {
      if (i == j) next

      d_ij <- distance_df$distance_meters[
        distance_df$county1 == i & distance_df$county2 == j
      ]
      if (length(d_ij) == 0 || is.na(d_ij)) { beta_mat[i, j] <- NA; next }

      # Counties k != i,j that are closer to i than j
      closer <- distance_df$county1 == i &
        distance_df$county2 != i &
        distance_df$county2 != j &
        distance_df$distance_meters < d_ij

      s_ij <- sum(N[distance_df$county2[closer]], na.rm = TRUE)
      Ni   <- N[i]; Nj <- N[j]

      beta_mat[i, j] <- (Ni * Nj) / ((Ni + s_ij) * (Ni + Nj + s_ij))
    }
  }
  return(beta_mat)
}


#' Run a single stochastic SIR simulation
#'
#' @param N_vec           Named numeric vector of patch populations (in patch order)
#' @param beta_between    Between-patch transmission rate matrix (n x n)
#' @param beta_within     Within-patch transmission rate (scalar)
#' @param gamma           Recovery rate
#' @param seed_patches    Integer vector of patch indices to seed with 1 infection each
#' @param time_steps      Number of discrete time steps
#' @return Data frame with columns: time, patch, S, I, R

run_sir_simulation <- function(N_vec, beta_between, beta_within,
                               gamma, seed_patches, time_steps) {
  n <- length(N_vec)

  # Initial conditions: seed 1 infected individual per seed patch
  I <- rep(0L, n)
  I[seed_patches] <- 1L
  S <- as.integer(N_vec) - I
  R <- rep(0L, n)

  results <- vector("list", time_steps)

  for (t in seq_len(time_steps)) {
    # Within-patch transmission
    new_inf <- rbinom(n, S, beta_within * I / N_vec)

    # Between-patch transmission
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) {
          new_inf[i] <- new_inf[i] +
            rbinom(1, S[i], beta_between[i, j] * I[j] / N_vec[i])
        }
      }
    }

    recoveries <- rbinom(n, I, gamma)

    # Update compartments
    S <- S - new_inf
    I <- I + new_inf - recoveries
    R <- R + recoveries

    results[[t]] <- data.frame(time = t, patch = seq_len(n), S = S, I = I, R = R)
  }

  return(do.call(rbind, results))
}


#' Compute mean first-case arrival time per patch across simulations
#'
#' @param sim_results   Combined simulation data frame (needs: simulation, patch, I, time)
#' @param patch_lookup  Data frame mapping patch index to county name
#' @return Data frame with patch, avg_first_case_time, county

compute_first_case_times <- function(sim_results, patch_lookup) {
  sim_results %>%
    filter(I > 0) %>%
    group_by(simulation, patch) %>%
    summarise(first_case_time = min(time), .groups = "drop") %>%
    group_by(patch) %>%
    summarise(avg_first_case_time = mean(first_case_time), .groups = "drop") %>%
    left_join(patch_lookup, by = "patch")
}


#' Plot effective distance vs. first-case arrival time
#'
#' @param df            Data frame with columns: distance, first_week_county2, county2
#' @param source_county Name of the source county (for title)
#' @param x_label       X-axis label
#' @return ggplot object

plot_distance_vs_onset <- function(df, source_county, x_label = "Effective distance") {
  ggplot(df, aes(x = distance, y = first_week_county2)) +
    geom_point(color = "#2c7fb8", size = 2.5) +
    geom_smooth(method = "loess", color = "#045a8d") +
    geom_text_repel(aes(label = county2), size = 4, max.overlaps = 36, color = "black") +
    labs(
      title = paste("Source:", source_county),
      x     = x_label,
      y     = "Average time for first case (time steps)"
    ) +
    theme_classic() +
    theme(
      axis.line  = element_line(linewidth = 1),
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5)
    )
}


###############################################################################
#                                                                             #
#                    MODEL 1: RADIATION MODEL APPROACH                        #
#                                                                             #
###############################################################################

cat("\n===================================================================\n")
cat("  MODEL 1: Radiation Model\n")
cat("===================================================================\n\n")

# ---- 6a. Parameters specific to radiation model ----
time_steps_rad    <- 100   # Simulation length
beta_within_rad   <- 0.3   # Within-patch transmission rate
seed_patches_rad  <- c(10, 23)  # Patch indices to seed (verify mapping below)

# Patch-to-county mapping (socio$county order)
patch_lookup_rad <- data.frame(patch = 1:n_patches, county = socio$county)

cat("Patch-to-county mapping:\n")
print(patch_lookup_rad, row.names = FALSE)
cat("\nSeed patches:\n")
for (sp in seed_patches_rad) {
  cat(sprintf("  Patch %d -> %s (pop: %s)\n", sp, socio$county[sp],
              format(N[sp], big.mark = ",")))
}

# ---- 6b. Compute radiation-based coupling matrix ----
cat("\nComputing radiation model coupling matrix...\n")
beta_rad_raw <- compute_radiation_classical(N, distance_df)

# Scale to a realistic between-patch transmission rate
beta_rad <- beta_rad_raw / max(beta_rad_raw, na.rm = TRUE) * 0.01

# Effective distance: d_eff = 1 - log(coupling)
# Used for correlation analysis (higher coupling -> shorter effective distance)
beta_rad_named <- beta_rad
colnames(beta_rad_named) <- rownames(beta_rad_named) <- socio$county

eff_dist_rad <- as.data.frame(as.table(beta_rad_named))
colnames(eff_dist_rad) <- c("county1", "county2", "coupling")
eff_dist_rad$distance <- 1 - log(eff_dist_rad$coupling)

# ---- 6c. Run simulations ----
cat(sprintf("Running %d simulations (%d time steps)...\n", n_simulations, time_steps_rad))
plan(multisession, workers = parallel::detectCores() - 2)

sim_results_rad <- future_lapply(1:n_simulations, function(sim) {
  df <- run_sir_simulation(
    N_vec        = N,
    beta_between = beta_rad,
    beta_within  = beta_within_rad,
    gamma        = gamma,
    seed_patches = seed_patches_rad,
    time_steps   = time_steps_rad
  )
  df$simulation <- sim
  return(df)
})

sim_results_rad <- do.call(rbind, sim_results_rad)

# ---- 6d. Compute first-case arrival times ----
first_case_rad <- compute_first_case_times(sim_results_rad, patch_lookup_rad)

cat("\nFirst case arrival times (radiation model):\n")
print(first_case_rad %>% arrange(avg_first_case_time), n = n_patches)

# ---- 6e. Merge distances with arrival times and plot ----
eff_dist_rad <- eff_dist_rad %>%
  left_join(first_case_rad, by = c("county2" = "county")) %>%
  rename(first_week_county2 = avg_first_case_time)

# Filter for a single source county (Greenville)
source_rad <- "Greenville"
df_plot_rad <- eff_dist_rad %>%
  filter(county1 == source_rad, county2 != source_rad)

fig_rad <- plot_distance_vs_onset(df_plot_rad, source_rad,
                                  x_label = "Effective distance (radiation model)")
print(fig_rad)

# Spearman correlation
cat("\nSpearman correlation (radiation model, source =", source_rad, "):\n")
print(cor.test(df_plot_rad$distance, df_plot_rad$first_week_county2, method = "spearman"))


###############################################################################
#                                                                             #
#                    MODEL 2: COMMUTING DATA APPROACH                         #
#                                                                             #
###############################################################################

cat("\n===================================================================\n")
cat("  MODEL 2: Commuting Data\n")
cat("===================================================================\n\n")

# ---- 7a. Load and prepare commuting data ----
commute_data <- read.csv("commute_SC_clean.csv")
commute_data <- commute_data %>%
  select(-X, -State_origin, -State_work)
commute_data$County_origin <- gsub(" County$", "", commute_data$County_origin)
commute_data$County_work   <- gsub(" County$", "", commute_data$County_work)

# Build wide-format adjacency matrix of commuting flows
adj_matrix <- commute_data %>%
  pivot_wider(
    names_from   = County_work,
    values_from  = Flow_volume,
    values_fill  = list(Flow_volume = 0)
  ) %>%
  as.data.frame()

rownames(adj_matrix) <- adj_matrix$County_origin
adj_matrix <- adj_matrix[, -1]

# Ensure consistent row/column ordering
common_order <- intersect(rownames(adj_matrix), colnames(adj_matrix))
adj_matrix   <- as.matrix(adj_matrix[common_order, common_order])
diag(adj_matrix) <- 0

# ---- 7b. Compute effective distance from commuting flows ----
# Transition probability: P_ij = flow_ij / total_outflow_i
P_ij <- adj_matrix / rowSums(adj_matrix, na.rm = TRUE)

# Effective distance via shortest path on -log(P_ij) graph
g_commute <- graph_from_adjacency_matrix(P_ij, mode = "directed",
                                          weighted = TRUE, diag = FALSE)
E(g_commute)$weight <- -log(E(g_commute)$weight)
P_eff <- distances(g_commute, mode = "out", weights = E(g_commute)$weight)

# Between-patch transmission matrix (scaled commuting probability)
beta_commute <- P_ij * 0.001

# ---- 7c. Patch ordering and parameters ----
# IMPORTANT: Patches in the commuting model follow common_order, which may
# differ from socio$county order. All vectors must be aligned accordingly.

patch_lookup_commute <- data.frame(patch = 1:length(common_order), county = common_order)
N_commute <- N[common_order]  # Population reordered to match commuting matrix

time_steps_commute   <- 120
beta_within_commute  <- 0.4
seed_patches_commute <- c(4, 35)

cat("Patch-to-county mapping (commuting order):\n")
print(patch_lookup_commute, row.names = FALSE)
cat("\nSeed patches:\n")
for (sp in seed_patches_commute) {
  cat(sprintf("  Patch %d -> %s (pop: %s)\n", sp, common_order[sp],
              format(N_commute[sp], big.mark = ",")))
}

# ---- 7d. Run simulations ----
cat(sprintf("\nRunning %d simulations (%d time steps)...\n",
            n_simulations, time_steps_commute))

sim_results_commute <- future_lapply(1:n_simulations, function(sim) {
  df <- run_sir_simulation(
    N_vec        = N_commute,
    beta_between = beta_commute,
    beta_within  = beta_within_commute,
    gamma        = gamma,
    seed_patches = seed_patches_commute,
    time_steps   = time_steps_commute
  )
  df$simulation <- sim
  return(df)
})

sim_results_commute <- do.call(rbind, sim_results_commute)

# ---- 7e. Compute first-case arrival times ----
first_case_commute <- compute_first_case_times(sim_results_commute, patch_lookup_commute)

cat("\nFirst case arrival times (commuting model):\n")
print(first_case_commute %>% arrange(avg_first_case_time), n = n_patches)

# ---- 7f. Merge distances with arrival times and plot ----
eff_dist_commute <- as.data.frame(as.table(P_eff))
colnames(eff_dist_commute) <- c("county1", "county2", "distance")

eff_dist_commute <- eff_dist_commute %>%
  left_join(first_case_commute, by = c("county2" = "county")) %>%
  rename(first_week_county2 = avg_first_case_time)

# Filter for a single source county (Greenville)
source_commute <- "Greenville"
df_plot_commute <- eff_dist_commute %>%
  filter(county1 == source_commute, county2 != source_commute)

fig_commute <- plot_distance_vs_onset(df_plot_commute, source_commute,
                                      x_label = "Effective distance (commuting)")
print(fig_commute)

# Spearman correlation
cat("\nSpearman correlation (commuting model, source =", source_commute, "):\n")
print(cor.test(df_plot_commute$distance, df_plot_commute$first_week_county2, method = "spearman"))

# =============================================================================
# 8. SESSION INFO (for reproducibility)
# =============================================================================
cat("\n===================================================================\n")
cat("  Session Info\n")
cat("===================================================================\n\n")
sessionInfo()
