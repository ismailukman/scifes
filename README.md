# SCIFES Graph Theory Experiment

## Overview

This repository contains the analysis pipeline for investigating dynamic functional brain network organization in Spinal Cord Injury (SCI) patients receiving Functional Electrical Stimulation (FES) therapy. Analysis based on collected rs-fRMI data, and graph based  multilayer community detection method..

## Experiment Pipeline

### Phase 1: Time-Resolved Functional Connectivity Extraction (Python)

**Script**: `SCIFES_Multilayer_by_TPs.ipynb`

#### Input Data

- **Format**: Raw resting-state fMRI NIfTI files (.nii/.nii.gz)
- **Groups**:
  - PreFES: 7 subjects (14 runs, 2 per subject)
  - PreNFES: 5 subjects (10 runs, 2 per subject)
  - PostFES: 7 subjects (14 runs, 2 per subject)
  - PostNFES: 5 subjects (10 runs, 2 per subject)
- **Imaging Parameters**: TR = 2s, ~420 TRs per subject (2 concatenated runs)

#### Processing Steps

1. **Regional Time Series Extraction**

   - Atlas: Schaefer 2018 (200 parcels)
   - Method: NiftiLabelsMasker extracts BOLD signals per region
   - Output: (time_points, 200) per run
2. **Run Concatenation**

   - Concatenates 2 runs per subject along time axis
   - Output: (~420 TRs, 200 regions) per subject
3. **Sliding Window Segmentation**

   - Window size: 30 seconds (15 TRs)
   - Step size: 8 TRs (configurable)
   - Windows per subject: 50
   - Output: List of windowed time series (15 TRs × 200 regions)
4. **Functional Connectivity Computation**

   - Method: Pearson correlation per window
   - Output per subject: (50 windows, 200 regions, 200 regions)

#### Output Files

Saved to `sci_data/SCI/fc/`:

- `corr_pre_fes_8tr_windows.mat` - Shape: (7, 50, 200, 200)
- `corr_pre_nfes_8tr_windows.mat` - Shape: (5, 50, 200, 200)
- `corr_post_fes_8tr_windows.mat` - Shape: (7, 50, 200, 200)
- `corr_post_nfes_8tr_windows.mat` - Shape: (5, 50, 200, 200)

Format: (n_subjects, n_windows, n_regions, n_regions)

---

### Phase 2: Multilayer Community Detection (MATLAB)

**Script**: `cons_mlcd_win_parallel_v3.m`
**Core Function**: `multilayer_community_detection_individual.m`

#### Input Data

- Loads FC tensors from Phase 1 (.mat files)
- Configuration: gamma=1.0, omega=1.0

#### Processing Steps

1. **Data Loading & Preprocessing**

   - Parallel loading of all 4 group FC tensors
   - Permutation to (200, 200, subjects×windows) format
   - Matrix symmetrization and cleaning (NaN/Inf handling)
   - Random window order shuffling per subject
2. **Multilayer Reshaping**

   - Reshapes from flat (subjects×windows) to subject-specific structure
   - Creates cell arrays: (n_subjects × n_windows_per_subject)
   - Each cell contains one (200×200) FC matrix
3. **Multilayer Community Detection**

   - **Method**: Ordinal multilayer modularity (`multiord`)
   - **Coupling**: Temporal coupling between consecutive windows (omega=1.0)
   - **Resolution**: Intra-layer resolution parameter (gamma=1.0)
   - **Iterations**: 100 repetitions per subject for consensus
   - **Algorithm**: `iterated_genlouvain` with `postprocess_ordinal_multilayer`
4. **Consensus Community Detection**

   - For each subject and window:
     - Computes **allegiance matrix** from 100 iterations
     - Detects communities on allegiance matrix
     - Iterates until all 100 assignments converge (consensus reached)
     - Matches community labels across windows for consistency
   - Outputs:
     - `multi_module_consensus`: Community assignments (200 nodes × 50 windows)
     - `multi_comm_consensus`: Number of communities per window
     - `multi_modQ`: Modularity quality scores
5. **Results Extraction**

   - Extracts consensus partitions for all subjects
   - Computes mode modularity values across iterations
   - Collects community counts per window/subject

#### Output Files

Saved to `sci_data/SCI/modularity_var/`:

- `mlcd_prefes_8tr_wins.mat`:

  - `N_all_g_prefes`: Consensus community assignments (200 × 350) [7 subjects × 50 windows]
  - `Q_g_prefes`: Modularity scores per subject
  - `comm_cons_all_g_prefes`: Number of communities per subject
- `mlcd_prenfes_8tr_wins.mat`: Same structure for PreNFES (5 subjects)
- `mlcd_postfes_8tr_wins.mat`: Same structure for PostFES (7 subjects)
- `mlcd_postnfes_8tr_wins.mat`: Same structure for PostNFES (5 subjects)

---

### Multilayer Community Detection

**Purpose**: Detect consensus community structure across temporal layers (windows)

**Algorithm**:

1. Constructs multilayer modularity matrix B using `multiord` (temporal coupling)
2. Runs Louvain algorithm 100 times via `iterated_genlouvain`
3. For each window:

   - Computes allegiance matrix (co-assignment frequency)
   - Detects communities on allegiance matrix
   - Iterates until consensus (all 100 iterations agree)

**Key Parameters**:

- `gamma`: Intra-layer resolution (controls community size within windows)
- `omega`: Inter-layer coupling (temporal coupling strength between windows)
- `coupling_type`: 'ord' for ordinal (temporal) or 'cat' for categorical
- `n_repeat`: 100 iterations for consensus stability

---

Copyright (C) 2025
