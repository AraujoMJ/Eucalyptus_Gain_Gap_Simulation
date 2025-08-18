# Eucalyptus_Gain_Gap_Simulation
R scripts and data for simulating genetic/genotypic gain gaps in Eucalyptus pellita under varying selfing rates. Includes pedigree, genomic, and single-step models with dominance, using the main function Tree_breeding_simulator.R to generate, analyze, and visualize results.

# Quantifying Genetic and Genotypic Gain Gaps in *Eucalyptus*

This repository contains simulation scripts, data, and documentation for the study:

**Araujo MJ, Bush D, Tambarussi EV. _Quantifying genetic and genotypic gain gaps in *Eucalyptus*: The Hidden Cost of Ignoring Inbreeding and Dominance_. Heredity.**

## Overview

The repository provides the code and resources needed to reproduce simulations and analyses quantifying the impact of inbreeding and dominance on genetic and genotypic gains in *Eucalyptus pellita* breeding programs.

The **main simulation pipeline** is implemented in:

- **`Tree_breeding_simulator.R`** – Generates simulated datasets, fits multiple models (pedigree, genomic, and single-step with dominance), and outputs genetic parameter estimates, selection gains, and gain gap statistics.

---

## Simulation Highlights

- **Species:** *Eucalyptus pellita*
- **Mating system:** Mixed, with controlled selfing rates (0–95%)
- **Models tested:**  
  - Pedigree-based (Half-Sib, Family models)  
  - Genomic (GBLUP, GDBLUP)  
  - Single-step genomic with dominance (ssGDBLUP)  
- **Effects modeled:** Additive and dominance variance
- **Deployment strategies:**  
  - Genetic gain (BV-based) for breeding  
  - Genotypic gain (BV+dominance) for clonal deployment
- **Outputs:** Genetic parameters, heritabilities, dominance coefficients, selection accuracies, gain gaps, tables, and figures

---

## Repository Structure

```
Eucalyptus_Gain_Gap_Simulation/
│
├── R/
│   ├── Tree_breeding_simulator.R   # Main custom simulation function 
│   ├── Genetic_Gain.R              # Function to calculate genetic/genotypic gains
│   ├── Selfing_Eucalyptus_2025.R   # Script applying Tree_breeding_simulator to run simulations
│
├── data/
│   ├── output/                     # Simulation results (.qs files)
│
├── README.md
└── LICENSE

```
---

## Requirements

- **R version:** ≥ 4.2  
- **Required packages:**  
  - `AlphaSimR`
  - `ASReml-R`
  - `tidyverse`
  - `qs`
  - `data.table`
  - `ggplot2`

---

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/USERNAME/Eucalyptus_Gain_Gap_Simulation.git
   cd Eucalyptus_Gain_Gap_Simulation
   ```
2. Open and edit simulation parameters in `Selfing_Eucalyptus_2025.R` (e.g., selfing rates, number of replicates, selection intensities).

3. Run the main pipeline in R:

```{r}
source("R/Tree_breeding_simulator.R")
source("R/Selfing_Eucalyptus_2025.R")

```

results <- Tree_breeding_simulator(params)

4. Generated output tables and figures will be stored in `data/output/` and `figs/`.

## Citation
If you use this repository, please cite:

[Araujo MJ](https://orcid.org/0000-0001-7218-3879), Bush D, [Tambarussi EV](https://orcid.org/0000-0001-9478-5379). Quantifying genetic and genotypic gain gaps in *Eucalyptus*: The Hidden Cost of Ignoring Inbreeding and Dominance. Heredity. [https://doi.org/10.1038/s41437-025-00792-8](https://doi.org/10.1038/s41437-025-00792-8)
