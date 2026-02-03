# Metapopulation SIR Model for South Carolina Counties

A stochastic SIR (Susceptible-Infected-Recovered) metapopulation model simulating disease spread across 46 South Carolina counties using two mobility frameworks.

## Overview

This repository implements a spatial disease transmission model that couples local (within-county) transmission dynamics with between-county mobility. Two approaches for estimating inter-county movement are compared:

1. **Radiation Model**: Estimates mobility fluxes based on population sizes and geographic distance, following the classical radiation model framework (Simini et al., 2012).

2. **Commuting Data Model**: Uses empirical county-to-county commuting flows to parameterize between-patch transmission.

Both models analyze the relationship between **effective distance** and **disease arrival time**, demonstrating how mobility patterns influence epidemic spread.

## Model Description

### SIR Dynamics

Each county is modeled as a patch with local SIR dynamics:
- **Within-patch transmission**: Frequency-dependent transmission with rate β = 0.4
- **Recovery rate**: γ = 0.25 (4-day average infectious period)
- **Between-patch transmission**: Scaled by mobility-derived coupling matrix

### Radiation Model

The radiation model calculates movement probability between counties i and j as:

```
T_ij = (N_i × N_j) / ((N_i + s_ij) × (N_i + N_j + s_ij))
```

Where:
- N_i, N_j = populations of origin and destination counties
- s_ij = total population within radius d_ij (excluding i and j)

### Effective Distance

For the commuting model, effective distance is computed as:

```
d_eff = -log(P_ij)
```

Where P_ij is the transition probability (normalized commute flow). Shortest paths through this effective distance network predict disease arrival times.

## Required Data Files

| File | Description |
|------|-------------|
| `SC_county_sociodemographic.csv` | County populations (columns: County, Total_population) |
| `sc_shape_files.RData` | County shapefiles as sf object (sc_county_sf) |
| `commute_SC_clean.csv` | County-to-county commute flows (columns: County_origin, County_work, Flow_volume) |

## Dependencies

```r
install.packages(c(
  "dplyr",
  "ggplot2", 
  "sf",
  "tidyr",
  "igraph",
  "future.apply",
  "ggrepel"
))
```

## Usage

1. Place required data files in your working directory
2. Update file paths in the script if necessary
3. Run the script:

```r
source("metapopulation_sir_model.R")
```

The script will:
- Run 1,000 stochastic simulations for each mobility model
- Calculate average first case arrival times for each county
- Generate scatter plots of effective distance vs. arrival time
- Perform Spearman correlation analysis

## Output

- **Figures**: Scatter plots showing the relationship between effective distance from the seed location (Greenville) and average time to first case for all other counties
- **Statistical analysis**: Spearman correlation coefficient and p-value
- **Optional CSV exports**: Full simulation results (uncomment `write.csv` lines)

## Customization

### Change seed location
Modify the `initial_infected` vector in the simulation functions to seed infections in different counties.

### Adjust parameters
Key parameters at the top of each section:
- `n_simulations`: Number of Monte Carlo runs (default: 1000)
- `time_steps`: Simulation duration (default: 100-120)
- `gamma`: Recovery rate (default: 0.25)
- Within-patch transmission rate (default: 0.4)
- Between-patch scaling factor (default: 0.001)

## References

- Simini, F., González, M. C., Maritan, A., & Barabási, A. L. (2012). A universal model for mobility and migration patterns. *Nature*, 484(7392), 96-100.
- Brockmann, D., & Helbing, D. (2013). The hidden geometry of complex, network-driven contagion phenomena. *Science*, 342(6164), 1337-1342.

## License

[Add your license here]

## Contact

[Add contact information here]
