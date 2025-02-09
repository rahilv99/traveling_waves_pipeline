# Attention Waves Analysis Pipeline

This repository contains a MATLAB-based pipeline for analyzing attention-related neural oscillations from intracranial EEG data.  The pipeline processes raw EEG data, identifies oscillatory clusters, performs circular-linear regression analysis, and generates visualizations of the results.

## Features

* **Data Organization:** Organizes raw EEG data into structured formats suitable for analysis.
* **Power Spectrum Analysis:** Computes power spectral densities using Morlet wavelets.  Handles multiple event types (presentation, recognition, free recall).
* **Oscillatory Cluster Identification:** Identifies spatially clustered electrodes exhibiting significant oscillatory activity within specified frequency bands.
* **Circular-Linear Regression:** Fits a circular-linear regression model to investigate the relationship between neural phase and spatial location.
* **Biowulf Swarm Integration:** Leverages Biowulf for parallel processing of computationally intensive tasks.
* **Data Visualization:** Generates visualizations of oscillatory activity and directional trends on brain surface plots.
* **Statistical Analysis:** Performs statistical comparisons between cued and uncued conditions using multiple comparisons correction.

## Usage

The pipeline consists of several MATLAB scripts that should be run sequentially:

1. `attn_sessions.m`: Creates a session structure containing behavioral information.
2. `attn_organize.m`: Loads and organizes raw EEG data for each session.
3. `biowolf_swarm_pow.m`: Generates a Biowulf swarm script for parallel power spectrum analysis using `attn_pow_biowolf_v1.m`.
4. `merge_power.m`: Merges the power spectrum results from Biowulf.
5. `find_clusters.m`: Identifies oscillatory clusters based on power spectral density.
6. `preprocessing.m`: Preprocesses data for circular-linear regression.  Generates a Biowulf swarm script for parallel processing using `attn_fitting_biowolf_v1.m`.
7. `merge_processed.m`: Merges the processed data from Biowulf.
8. `wave_extraction_v2.m`: Extracts and analyzes wave information, generating summary structures.
9. `location_filter_1.m` and `location_filter_2.m`: Filter data based on location (ATL, FL, PL, PTL) and generate plots.
10. `avg_brain_plotter.m`: Creates average brain plots for visualization.
11. `beh_filters.m`: Filters patients based on behavioral performance criteria.

## Installation

1. Ensure MATLAB is installed.
2. Clone this repository.
3. Add necessary toolboxes and functions.  See "Dependencies".
4. Configure paths to data and output directories. See "Configuration".


## Technologies Used

* **MATLAB:**  The primary programming language for the entire pipeline.
* **Biowulf:** High-performance computing cluster for parallel processing.
* **Morlet Wavelets:** Used for time-frequency analysis of EEG data.
* **Circular Statistics Toolbox:** Used for circular-linear regression analysis.
* **PCA:** Principal Component Analysis is used to determine the electrode plane for circular-linear regression.


## Statistical Analysis

The pipeline employs multiple statistical methods:

* **Circular mean:** Used to compute the average phase angle.
* **Circular-linear regression:** Used to model the relationship between phase and spatial location.
* **Multiple Comparisons Correction:**  Used to control for false positives in statistical comparisons between cued and uncued conditions (using `stat_multcomp` function).

## Configuration

Modify the following parameters in the respective scripts:

* File paths for EEG data, behavioral data, and output directories.
* Parameters for power spectral density calculation (e.g., frequency range, wavelet width).
* Parameters for cluster identification (e.g., distance threshold, minimum number of electrodes).
* Parameters for circular-linear regression (e.g., distance threshold).

## Dependencies

* MATLAB Signal Processing Toolbox
* MATLAB Statistics and Machine Learning Toolbox
* Circular Statistics Toolbox for Matlab (circstat) - [Link to toolbox if available]
* Custom functions located in the `/Volumes/Rahil_FRNU/Scripts/ZaghloulCodebase` and `/Volumes/Rahil_FRNU/Scripts/(7) Analysis/functions` directories. You may request the owner of the repository for access.
