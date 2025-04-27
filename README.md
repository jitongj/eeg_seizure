# EEG Seizure Analysis and Tensor-Based Feature Extraction

This repository provides a comprehensive MATLAB-based framework for EEG seizure detection, prediction, and tensor decomposition analysis, utilizing the [CHB-MIT Scalp EEG Database](https://physionet.org/content/chbmit/1.0.0/).

It covers two main pipelines:

- Seizure Feature Extraction and Tensor Construction for classification tasks (detection and prediction)
- Seizure Tensor Decomposition Analysis for understanding spatial, temporal, and frequency patterns

## Repository Structure and Pipelines

| Task | Folder / Script | Method |
|:---|:---|:---|
| Basic EEG Tensor Generation | BCI/get_tensor.m | Simple window-based tensor extraction |
| Seizure Detection and Prediction (CWT) | CHB/build_all_tensors.m with chb_fun/ | Continuous Wavelet Transform (CWT) |
| Seizure Detection and Prediction (FFT) | CHB-FFT/CHB_fft.m with function/ | Fast Fourier Transform (FFT) |
| Seizure Detection and Prediction (Feature-based) | Feature_tensors/build_feature_tensor.m with function/ | Statistical features (mean, variance, energy) |
| Flexible Experiments (Custom Tensor Generation) | flexible_fun/build_all_tensors.m | Customized tensor generation pipelines |
| Tensor Decomposition Analysis | tensor_analysis/ (chb01.m, chb02.m, chb03.m) | CANDECOMP/PARAFAC (CPD) tensor decomposition |

## Folder Descriptions

### 1. BCI/
Contains simple tensor generation scripts for BCI Competition III EEG data.

- Main script: get_tensor.m
- Focus: Basic feature tensor construction for model training.

### 2. CHB/
The core folder for time-frequency tensor generation using Continuous Wavelet Transform (CWT) on CHB-MIT EEG data.

- Main script: build_all_tensors.m
- Subfolder: chb_fun/
  - EEG reading (read.m, read_eeg.m, read_seizure.m, read_alldata.m, read_preictal.m)
  - Label generation (label_seizure.m, label_alldata.m, label_eeg.m)
  - Filtering (filter_dc.m, filter_60Hz.m)
  - Window slicing (cut_window.m, cut_window_training.m, cut_window_validation.m)
  - Tensor construction (seizure_prediction_tensor.m, nonseizure_prediction_tensor.m)
  - Utility tools (clear_empty.m, perfcurve.m)
  - Channel adjustment (delete_chan.m)

### 3. CHB-FFT/
Focused on FFT-based time-frequency tensor generation.

- Main script: CHB_fft.m
- Subfolder: function/
  - EEG reading and labeling
  - Filtering
  - FFT-based tensor construction (fft_SeiPred_tensor.m, fft_SeiDet_tensor.m, fft_NonseiPred_tensor.m, fft_NonseiDet_tensor.m)

Note: Functions like filter_dc.m, label_seizure.m are duplicated here to make the FFT module fully self-contained.

### 4. Feature_tensors/
Focused on statistical feature extraction (mean, variance, energy) from EEG windows for seizure detection and prediction.

- Main script: build_feature_tensor.m
- Subfolder: function/
  - Feature extraction scripts (det_feature.m, pre_feature.m)
  - EEG reading and labeling
  - Filtering
  - Window slicing utilities

Note: Basic utility functions are duplicated locally for modularity.

### 5. flexible_fun/
Provides a flexible framework for customized tensor generation.

- Main script: build_all_tensors.m (formerly build_all_tensors_1.m, renamed for consistency)
- Flexible scripts:
  - s_pre.m, s_det.m (seizure prediction and detection tensors)
  - ns_pred.m, ns_det.m (nonseizure prediction and detection tensors)

This folder is useful for exploratory experiments that deviate from standard pipelines.

### 6. function/
A standalone toolbox collecting general-purpose helper functions.

- EEG reading and labeling
- Filtering (filter_dc.m, filter_60Hz.m)
- Channel deletion (delete_channels.m, delete_channels_1.m)
- Window slicing
- Performance evaluation (perfcurve.m)
- Miscellaneous utilities (statset.m)

Note: Although similar utilities exist in each module's local function/ folders, this folder provides a centralized version of key functions for broader usage.

### 7. tensor_analysis/
Contains tensor decomposition analysis for selected CHB-MIT samples.

- Scripts:
  - chb01.m
  - chb02.m
  - chb03.m
- Process:
  - Load .set EEG files
  - Preprocess and apply Continuous Wavelet Transform (CWT)
  - Conduct CPD tensor decomposition
  - Visualize spatial, temporal, and frequency patterns

Input Data:
Selected .set files: chb01_03.set, chb02_16+.set, chb03_01.set

## About Redundancy

Many basic helper functions like EEG reading, filtering, and labeling exist in each module's local function/ folder and in the centralized top-level function/ folder.

This duplication was intentional to:

- Ensure each module is fully self-contained and runnable independently
- Avoid complicated cross-folder dependencies

No urgent need for deduplication unless a major refactor is planned.

Typical duplicated functions include:

- filter_dc.m
- filter_60Hz.m
- label_seizure.m
- label_alldata.m
- read.m, read_eeg.m, read_seizure.m, etc.

## Data Processing and Outputs

- Seizure Detection:
  - Scripts: eeg_pro.m, seizure_pro.m
  - Output: 3D tensors
  - Each patient has 20 windows.

- Seizure Prediction:
  - Scripts: pre_eeg.m, pre_seizure.m
  - Output: 4D tensors
  - Each seizure is preceded by a 3.5-minute prediction window.
  - Number of nonseizure prediction windows matches the number of seizures per patient.

Supported Patients: 16 patients (patients 1–11, 14, and 20–23)

## Input Data

- CHB-MIT .edf or .set EEG recordings.
- For tensor decomposition (tensor_analysis/):
  - chb01_03.set
  - chb02_16+.set
  - chb03_01.set

Important:
Text summaries for all patients previously used for seizure timing (summary/ folder) are not included. Access must be requested separately if needed.

## Requirements

- MATLAB
- EEGLAB toolbox (for .set files)
- Tensor Toolbox (for tensor decomposition)

## Quick Start Example

```matlab
% Example: Run the CWT-based tensor construction for patient 1
build_all_tensors(1)
```
