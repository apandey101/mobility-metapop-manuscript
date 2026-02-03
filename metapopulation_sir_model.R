# =============================================================================
# Metapopulation SIR Model for South Carolina Counties
# =============================================================================
# This script implements a stochastic SIR (Susceptible-Infected-Recovered) 
# metapopulation model using two different mobility frameworks:
#   1. Classical Radiation Model (based on geographic distance and population)
#   2. Commuting Data (based on actual county-to-county commute flows)
#
# The model simulates disease spread across 46 South Carolina counties and
# analyzes the relationship between effective distance and disease arrival time.
# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(sf)
library(tidyr)
library(igraph)
library(future.apply)
library(ggrepel)

# =============================================================================
# PART 1: RADIATION MODEL FRAMEWORK
# =============================================================================
# The radiation model estimates mobility fluxes between locations based on
# population sizes and intervening opportunities (populations within radius).
# Reference: Simini et al. (2012) Nature
# =============================================================================

rm(list = ls())

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
# Load sociodemographic data for SC counties
# Required columns: County, Total_population
socio <- read.csv("SC_county_sociodemographic.csv")

# Extract population vector and rename county column
N <- as.numeric(socio$Total_population)
socio <- socio %>%
  rename(county = County)
names(N) <- socio$county

# Load SC county shapefiles
# Contains: sc_county_sf (sf object with county geometries)
load("sc_shape_files.RData")

# -----------------------------------------------------------------------------
# Compute Geographic Distance Matrix
# -----------------------------------------------------------------------------
# Transform to projected coordinate system for accurate distance calculations
counties <- st_transform(sc_county_sf, crs = 3857)  # Web Mercator projection

# Compute county centroids
county_centroids <- st_centroid(counties)

# Compute pairwise distances between centroids (in meters)
dist_matrix <- st_distance(county_centroids)
colnames(dist_matrix) <- rownames(dist_matrix) <- counties$NAME

# Convert distance matrix to long format
distance_df <- as.data.frame(as.table(dist_matrix))
colnames(distance_df) <- c("county1", "county2", "distance_meters")
distance_df$distance_meters <- as.numeric(distance_df$distance_meters)

# -----------------------------------------------------------------------------
# Model Parameters
# -----------------------------------------------------------------------------
n_patches <- 46        # Number of SC counties
gamma <- 0.25          # Recovery rate (1/infectious period)
time_steps <- 100      # Simulation duration
n_simulations <- 1000  # Number of stochastic simulations

# -----------------------------------------------------------------------------
# Radiation Model Implementation
# -----------------------------------------------------------------------------
#' Compute transmission rates using the Classical Radiation Model
#' 
#' The radiation model calculates the probability of movement between locations
#' based on: T_ij = (N_i * N_j) / ((N_i + s_ij) * (N_i + N_j + s_ij))
#' where s_ij is the total population within distance d_ij (excluding i and j)
#'
#' @param N Named vector of population sizes
#' @param distance_df Data frame with columns: county1, county2, distance_meters
#' @return Matrix of transmission rates between patches

compute_radiation_classical <- function(N, distance_df) {
  counties <- names(N)
  n_patches <- length(counties)
  
  beta_matrix <- matrix(0, n_patches, n_patches)
  rownames(beta_matrix) <- colnames(beta_matrix) <- counties
  
  for (i in counties) {
    for (j in counties) {
      if (i != j) {
        # Get direct distance between counties i and j
        d_ij <- distance_df$distance_meters[
          distance_df$county1 == i & distance_df$county2 == j
        ]
        
        if (length(d_ij) == 0 || is.na(d_ij)) {
          beta_matrix[i, j] <- NA
          next
        }
        
        # Identify counties k within distance d_ij (excluding i and j)
        relevant_rows <- distance_df$county1 == i &
          distance_df$county2 != i &
          distance_df$county2 != j &
          distance_df$distance_meters < d_ij
        
        k_names <- distance_df$county2[relevant_rows]
        
        # Compute s_ij: total intervening population
        s_ij <- sum(N[k_names], na.rm = TRUE)
        
        # Classical radiation model formula
        Ni <- N[i]
        Nj <- N[j]
        beta_matrix[i, j] <- (Ni * Nj) / ((Ni + s_ij) * (Ni + Nj + s_ij))
      }
    }
  }
  
  return(beta_matrix)
}

# Compute and normalize between-patch transmission rates
beta_between_matrix <- compute_radiation_classical(N, distance_df)
beta_between_matrix <- beta_between_matrix / max(beta_between_matrix) * 0.001

# -----------------------------------------------------------------------------
# Stochastic SIR Simulation Function (Radiation Model)
# -----------------------------------------------------------------------------
#' Run a single stochastic SIR simulation
#' 
#' Uses binomial draws for infection and recovery events
#' Initial infections seeded in specific counties (indices 10, 24, 23)
#' 
#' @return Data frame with time series of S, I, R for each patch

run_sir_simulation_radiation <- function() {
  # Initialize compartments
  # Seed initial infections in select counties
  initial_infected <- c(rep(0, 9), 1, rep(0, 12), 1, rep(0, 23))
  S <- N - initial_infected
  I <- c(rep(0, 22), 1, rep(0, 23))
  R <- rep(0, 46)
  
  results <- data.frame()
  
  for (t in 1:time_steps) {
    # Within-patch transmission (frequency-dependent)
    new_infections <- rbinom(n_patches, S, 0.4 * I / N)
    
    # Recovery events
    recoveries <- rbinom(n_patches, I, gamma)
    
    # Between-patch transmission (radiation model-based)
    for (i in 1:n_patches) {
      for (j in 1:n_patches) {
        if (i != j) {
          imported_cases <- rbinom(1, S[i], beta_between_matrix[i, j] * I[j] / N[i])
          new_infections[i] <- new_infections[i] + imported_cases
        }
      }
    }
    
    # Update compartments
    S <- S - new_infections
    I <- I + new_infections - recoveries
    R <- R + recoveries
    
    # Store results
    results <- rbind(results, data.frame(
      time = t,
      patch = 1:n_patches,
      S = S,
      I = I,
      R = R
    ))
  }
  
  return(results)
}

# -----------------------------------------------------------------------------
# Run Parallel Simulations (Radiation Model)
# -----------------------------------------------------------------------------
plan(multisession, workers = parallel::detectCores() - 2)

simulation_results_radiation <- future_lapply(1:n_simulations, function(sim) {
  df <- run_sir_simulation_radiation()
  df$simulation <- sim
  return(df)
})

simulation_results_radiation <- do.call(rbind, simulation_results_radiation)

# Optional: Save results
# write.csv(simulation_results_radiation, "simulation_results_radiation.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Analyze First Case Arrival Times (Radiation Model)
# -----------------------------------------------------------------------------
# Calculate average time to first infection for each patch
first_case_times_radiation <- simulation_results_radiation %>%
  filter(I > 0) %>%
  group_by(simulation, patch) %>%
  summarise(first_case_time = min(time), .groups = "drop") %>%
  group_by(patch) %>%
  summarise(avg_first_case_time = mean(first_case_time), .groups = "drop")

first_case_times_radiation$county <- socio$county

# Prepare data for effective distance analysis
colnames(beta_between_matrix) <- rownames(beta_between_matrix) <- socio$county

distance_df_analysis <- as.data.frame(as.table(beta_between_matrix))
colnames(distance_df_analysis) <- c("county1", "county2", "beta")

distance_df_analysis <- distance_df_analysis %>%
  left_join(first_case_times_radiation, by = c("county2" = "county")) %>%
  rename(first_case_time = avg_first_case_time)

# Filter to analyze spread from Greenville
distance_df_greenville_radiation <- distance_df_analysis %>%
  filter(county1 == "Greenville" & county2 != "Greenville") %>%
  mutate(effective_distance = 1 - log(beta))  # Transform to effective distance

# -----------------------------------------------------------------------------
# Visualization: Effective Distance vs First Case Time (Radiation Model)
# -----------------------------------------------------------------------------
p1 <- ggplot(distance_df_greenville_radiation, 
             aes(x = effective_distance, y = first_case_time)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = 'lm', se = TRUE, color = "blue") +
  geom_text_repel(aes(label = county2), size = 3, max.overlaps = 40) +
  labs(
    title = "Radiation Model: Disease Spread from Greenville",
    x = "Effective Distance (1 - log(β))",
    y = "Average Time to First Case (days)"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  )

print(p1)

# =============================================================================
# PART 2: COMMUTING DATA FRAMEWORK
# =============================================================================
# Uses actual county-to-county commuting flows to estimate mobility
# and effective distances between counties.
# =============================================================================

rm(list = ls())

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
socio <- read.csv("SC_county_sociodemographic.csv")
N <- as.numeric(socio$Total_population)
socio <- socio %>%
  rename(county = County)
names(N) <- socio$county

# -----------------------------------------------------------------------------
# Model Parameters
# -----------------------------------------------------------------------------
n_patches <- 46
gamma <- 0.25
time_steps <- 120
n_simulations <- 1000

# -----------------------------------------------------------------------------
# Process Commuting Data
# -----------------------------------------------------------------------------
# Load commute data
# Required columns: County_origin, County_work, Flow_volume
commute_data <- read.csv("commute_SC_clean.csv")

# Clean county names (remove " County" suffix if present)
commute_data <- commute_data %>%
  select(-any_of(c("X", "State_origin", "State_work")))

commute_data$County_origin <- gsub(" County$", "", commute_data$County_origin)
commute_data$County_work <- gsub(" County$", "", commute_data$County_work)

# Create adjacency matrix from commute flows
adj_matrix <- commute_data %>%
  pivot_wider(
    names_from = County_work,
    values_from = Flow_volume,
    values_fill = list(Flow_volume = 0)
  ) %>%
  as.data.frame()

rownames(adj_matrix) <- adj_matrix$County_origin
adj_matrix <- adj_matrix[, -1]

# Ensure symmetric ordering of rows and columns
common_order <- intersect(rownames(adj_matrix), colnames(adj_matrix))
adj_matrix <- as.matrix(adj_matrix[common_order, common_order])

# Remove self-loops
diag(adj_matrix) <- 0

# -----------------------------------------------------------------------------
# Compute Effective Distance Matrix
# -----------------------------------------------------------------------------
# Transition probability matrix (row-normalized)
P_ij <- adj_matrix / rowSums(adj_matrix, na.rm = TRUE)

# Build directed weighted graph
g <- graph_from_adjacency_matrix(
  P_ij,
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)

# Transform weights to effective distance: d_eff = -log(P_ij)
E(g)$weight <- -log(E(g)$weight)

# Compute shortest path effective distances
P_eff <- distances(g, mode = "out", weights = E(g)$weight)

# Between-patch transmission rates (scaled transition probabilities)
beta_between_matrix <- P_ij * 0.001

# -----------------------------------------------------------------------------
# Stochastic SIR Simulation Function (Commute Data)
# -----------------------------------------------------------------------------
run_sir_simulation_commute <- function() {
  # Initialize compartments
  # Seed initial infection in county at index 4 (adjust as needed)
  initial_infected <- c(rep(0, 3), 1, rep(0, 42))
  S <- N - initial_infected
  I <- initial_infected
  R <- rep(0, 46)
  
  results <- data.frame()
  
  for (t in 1:time_steps) {
    # Within-patch transmission
    new_infections <- rbinom(n_patches, S, 0.4 * I / N)
    
    # Recovery events
    recoveries <- rbinom(n_patches, I, gamma)
    
    # Between-patch transmission (commute-based)
    for (i in 1:n_patches) {
      for (j in 1:n_patches) {
        if (i != j) {
          imported_cases <- rbinom(1, S[i], beta_between_matrix[i, j] * I[j] / N[i])
          new_infections[i] <- new_infections[i] + imported_cases
        }
      }
    }
    
    # Update compartments
    S <- S - new_infections
    I <- I + new_infections - recoveries
    R <- R + recoveries
    
    # Store results
    results <- rbind(results, data.frame(
      time = t,
      patch = 1:n_patches,
      S = S,
      I = I,
      R = R
    ))
  }
  
  return(results)
}

# -----------------------------------------------------------------------------
# Run Parallel Simulations (Commute Data)
# -----------------------------------------------------------------------------
plan(multisession, workers = parallel::detectCores() - 2)

simulation_results_commute <- future_lapply(1:n_simulations, function(sim) {
  df <- run_sir_simulation_commute()
  df$simulation <- sim
  return(df)
})

simulation_results_commute <- do.call(rbind, simulation_results_commute)

# Optional: Save results
# write.csv(simulation_results_commute, "simulation_results_commute.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Analyze First Case Arrival Times (Commute Data)
# -----------------------------------------------------------------------------
first_case_times_commute <- simulation_results_commute %>%
  filter(I > 0) %>%
  group_by(simulation, patch) %>%
  summarise(first_case_time = min(time), .groups = "drop") %>%
  group_by(patch) %>%
  summarise(avg_first_case_time = mean(first_case_time), .groups = "drop")

first_case_times_commute$county <- common_order

# Prepare effective distance analysis
distance_df_commute <- as.data.frame(as.table(P_eff))
colnames(distance_df_commute) <- c("county1", "county2", "effective_distance")

distance_df_commute <- distance_df_commute %>%
  left_join(first_case_times_commute, by = c("county2" = "county")) %>%
  rename(first_case_time = avg_first_case_time)

# Filter for Greenville as source
distance_df_greenville_commute <- distance_df_commute %>%
  filter(county1 == "Greenville" & county2 != "Greenville")

# -----------------------------------------------------------------------------
# Visualization: Effective Distance vs First Case Time (Commute Data)
# -----------------------------------------------------------------------------
p2 <- ggplot(distance_df_greenville_commute, 
             aes(x = effective_distance, y = first_case_time)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = 'lm', se = TRUE, color = "blue") +
  geom_text_repel(aes(label = county2), size = 3, max.overlaps = 40) +
  labs(
    title = "Commute Model: Disease Spread from Greenville",
    x = "Effective Distance (-log(P_ij) shortest path)",
    y = "Average Time to First Case (days)"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  )

print(p2)

# -----------------------------------------------------------------------------
# Statistical Analysis: Correlation Test
# -----------------------------------------------------------------------------
# Test correlation between effective distance and arrival time
cor_test <- cor.test(
  distance_df_greenville_commute$effective_distance,
  distance_df_greenville_commute$first_case_time,
  method = "spearman"
)

cat("\n=== Spearman Correlation Test (Commute Model) ===\n")
cat("Correlation coefficient (rho):", round(cor_test$estimate, 3), "\n")
cat("P-value:", format(cor_test$p.value, scientific = TRUE, digits = 3), "\n")

# =============================================================================
# END OF SCRIPT
# =============================================================================
