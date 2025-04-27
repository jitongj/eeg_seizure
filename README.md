# EEG Seizure Analysis Using Tensor Decomposition and Feature Extraction

This repository provides a complete framework for EEG seizure analysis using tensor decomposition, CWT, FFT, and statistical feature extraction methods. The data is primarily based on the [CHB-MIT Scalp EEG Database](https://physionet.org/content/chbmit/1.0.0/).

## Overview

The pipeline processes EEG data through multiple stages:
- Preprocessing (filtering, channel cleaning)
- Segmentation (window slicing around seizure and non-seizure events)
- Feature Extraction (Continuous Wavelet Transform, FFT, and statistical features)
- Tensor Construction for machine learning model training and analysis

The repository supports both:
- **Seizure detection** (identifying ongoing seizures)
- **Seizure prediction** (forecasting upcoming seizures)

The current implementation covers **16 patients**, specifically patients 1–11, 14, and 20–23.

**Prediction Tasks:**
- `pre_eeg.m` and `pre_seizure.m` generate 4D tensors for non-seizure and seizure prediction windows (3.5 minutes prior to each seizure).
- The number of prediction windows matches the number of seizures per patient.

**Detection Tasks:**
- `eeg_pro.m` and `seizure_pro.m` generate 3D tensors for non-seizure and seizure detection windows.
- Each patient has 20 windows (can be modified based on needs) for detection.

## Repository Structure

### Top-Level Folders

| Folder | Description |
|:---|:---|
| `BCI/` | Scripts for basic tensor generation from BCI Competition III EEG data. |
| `CHB/` | Core folder for time-frequency tensor generation using Continuous Wavelet Transform (CWT). |
| `CHB-FFT/` | Tensor generation using Fast Fourier Transform (FFT) features. |
| `feature_tensors/` | Tensor construction using statistical feature extraction from EEG windows. |
| `flexible_fun/` | Framework for customized flexible tensor generation (detection/prediction windows). |
| `function/` | Centralized toolbox collecting all basic helper functions (EEG reading, filtering, slicing, etc.). |
| `summary/` | Contains `summary.xlsx`, a summary of seizure start times per patient. (Note: original text summary files are not publicly included and require permission.) |
| `tensor_analysis/` | Contains `chb01.m`, `chb02.m`, `chb03.m` scripts analyzing selected EEG samples (`chb01_03.set`, `chb02_16+.set`, `chb03_01.set`) via CWT and tensor decomposition. |

### Folder Details

#### 1. `BCI/`
Simple tensor generation scripts for BCI EEG data.

- Main Script: `get_tensor.m`
- Focus: Basic feature tensor construction for model training.

#### 2. `CHB/`
The core folder for CWT-based time-frequency tensor generation using CHB-MIT EEG data.

- Main Script: `build_all_tensors.m`
- Subfolder: `chb_fun/`
  - EEG reading: `read.m`, `read_eeg.m`, `read_seizure.m`, `read_alldata.m`, `read_preictal.m`
  - Label generation: `label_seizure.m`, `label_alldata.m`, `label_eeg.m`
  - Filtering: `filter_dc.m`, `filter_60Hz.m`
  - Window slicing: `cut_window.m`, `cut_window_training.m`, `cut_window_validation.m`
  - Tensor creation: `seizure_prediction_tensor.m`, `nonseizure_prediction_tensor.m`
  - Utility tools: `clear_empty.m`, `perfcurve.m`
  - Channel adjustment: `delete_chan.m`

#### 3. `CHB-FFT/`
FFT-based tensor generation pipeline.

- Main Script: `CHB_fft.m`
- Subfolder: `function/`
  - EEG reading, labeling, filtering
  - FFT-based tensor construction: `fft_SeiPred_tensor.m`, `fft_SeiDet_tensor.m`, `fft_NonseiPred_tensor.m`, `fft_NonseiDet_tensor.m`

Note: Basic functions such as `filter_dc.m`, `label_seizure.m` are duplicated here for modularity.

#### 4. `feature_tensors/`
Statistical feature extraction from EEG windows for seizure detection and prediction.

- Main Script: `build_feature_tensor.m`
- Subfolder: `function/`
  - Feature extraction: `det_feature.m`, `pre_feature.m`
  - EEG reading, labeling, filtering
  - Window slicing utilities

Note: Contains duplicate helper functions to ensure folder independence.

#### 5. `flexible_fun/`
Customized tensor generation framework for flexible experimental settings.

- Main Script: `build_all_tensors.m` 
- Supporting scripts:
  - `s_pre.m`, `s_det.m`: Seizure prediction/detection tensor generation
  - `ns_pred.m`, `ns_det.m`: Nonseizure prediction/detection tensor generation

Useful for manually controlling window parameters or patient selections.

#### 6. `function/`
Standalone toolbox collecting important helper functions.

- EEG reading
- Label generation
- Filtering (`filter_dc.m`, `filter_60Hz.m`)
- Channel deletion (`delete_channels.m`, `delete_channels_1.m`)
- Window slicing
- Performance evaluation (`perfcurve.m`)
- Miscellaneous utilities (`statset.m`)

Note: This folder centralizes all necessary basic functions, but similar versions exist locally in other folders for self-contained usage.

#### 7. `summary/`
- File: `summary.xlsx`
- Description: Summarizes seizure start times across all selected patients.
- Note: The raw `.txt` summaries used for generating this file are not publicly available. Access requires permission.

#### 8. `tensor_analysis/`
Tensor decomposition analysis scripts for selected CHB-MIT EEG records.

- Files:
  - `chb01.m`: Analysis of `chb01_03.set`
  - `chb02.m`: Analysis of `chb02_16+.set`
  - `chb03.m`: Analysis of `chb03_01.set`
- Workflow:
  - Preprocessing with FIR filter
  - Tensor construction using CWT
  - CPD decomposition
  - Rank analysis and uniqueness verification
  - Visualization of spatial, temporal, and frequency patterns

### Redundancy Note

Many basic helper functions are intentionally duplicated across:
- `chb_fun/`
- `CHB-FFT/function/`
- `feature_tensors/function/`
- `function/`

This duplication is to keep each major module fully self-contained. No immediate refactoring is necessary unless a major reorganization is planned.

Typical duplicated files:
- `filter_dc.m`
- `filter_60Hz.m`
- `label_seizure.m`
- `label_alldata.m`
- `read.m`, `read_eeg.m`, `read_seizure.m`

---

## Requirements

- MATLAB
- EEGLAB Toolbox (for loading `.set` files)
- Tensor Toolbox (for CPD decomposition)

## Notes

- Tensor decomposition is repeated with multiple initializations to verify uniqueness.
- The selected rank for decomposition is based on consistency and error evaluation.

