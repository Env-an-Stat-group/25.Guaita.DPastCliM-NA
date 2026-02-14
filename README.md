# Downscaling Pipeline

This repository implements a statistical downscaling workflow using Principal Component Regression (PCR) (and optionally Linear Regression, LR). The pipeline consists of five main stages: data preparation, preprocessing, model calibration, downscaling and gridding, and visualization.

---

## 1. Data Preparation

### 1.1 Observational Data
Download GHCN-m observational data and store in:
<main_path>/obs_data

### 1.2 CMIP Data
Download CMIP model outputs and store in:
<main_path>/CMIP_data

### 1.3 Static Maps
Prepare static maps for the downscaling grid:
- Latitude array  
- Longitude array  
- Land–sea mask (IMERG is used)

Store in:
<main_path>/static_maps

---

## 2. Preprocessing

Run:
<main_path>/matlab_code/preprocessing/GHCNm.m

This script preprocesses the observational data.

---

## 3. Model Calibration

Run:
<main_path>/matlab_code/downscaling/PCR_calibration_v5.m

This calibrates the PCR model month by month.

Outputs are saved in:
<main_path>/downscaling_models_[model_name]

Files:
- Calibrated model:
  [variable]_PCR_models_mth=[month]_[domain_suffix]

- Original observation table:
  [variable]_PCR_original_data_mth=[month]_[domain_suffix]

---

## 4. Downscaling ESM Products

Directory:
<main_path>/matlab_code/downscaling

### 4.1 Station-Level Downscaling

Run:
PCR_downscaling_v5_stautocorr.m

(Optional Linear Regression):
PCR_downscaling_v5_stautocorr_LR.m

Outputs are saved in:
<main_path>/downscaling_output_[model_name]

File:
[name_var]_[model_name]_[experiment_name]_raw_PCRdownscaled_[domain_suffix]

This file contains reconstructed time series at each station.

### 4.1.1 Postprocessing of station-Level Downscaling

Run:
postprocessing/output_realizations.m


---

### 4.2 Gridding

Run:
PCR_downscaling_gridding_v1.m  
PCR_downscaling_gridding_v1_obs.m  

(Optional Linear Regression):
PCR_downscaling_gridding_v1_LR.m

Interpolation method: natural neighbor.

Outputs in <main_path>/downscaling_output_[model_name]:

i. Gridded downscaled product:
[variable]_[model_name]_[experiment_name]_PCRdownscaled_[domain_suffix]_[start_month_window]_[end_month_window]

ii. Prediction interval:
[variable]_[model_name]_[experiment_name]_PI_PCRdownscaled_[domain_suffix]_[start_month_window]_[end_month_window]

iii. Gridded observations:
[variable]_observations_[domain_suffix]

Contains:
- obsValue_map (all stations)  
- obsValue_test_map (testing stations only)  
- tgt_ESM (regridded ESM product)

iv. Metadata:
[variable]_metadata_table_[domain_suffix].csv

---

## 5. Visualization

Run:
<main_path>/matlab_code/visualization/PCR_plots_v5.m

This script generates maps, time series plots, and performance diagnostics.
