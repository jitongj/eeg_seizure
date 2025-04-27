# EEG Seizure Analysis Using Tensor Decomposition

This repository contains MATLAB scripts for analyzing EEG signals from the [CHB-MIT Scalp EEG Database](https://physionet.org/content/chbmit/1.0.0/).

We selected three samples:
- `chb01_03.set`
- `chb02_16+.set`
- `chb03_01.set`

## Overview

For each EEG record, we:
- Pre-processed the data using a basic FIR filter (band-pass: 0.5–40 Hz).
- Extracted 5-second and 10-second EEG tensors from both non-seizure and seizure periods.
- Standardized the signals per channel.
- Applied Continuous Wavelet Transform (CWT) to capture time–frequency features.
- Conducted CANDECOMP/PARAFAC (CPD) tensor decomposition.
- Evaluated the low-rank property and uniqueness of the tensor decompositions.
- Visualized the extracted spatial, temporal, and frequency patterns for seizure and non-seizure segments.

## Scripts

| Script | Description |
|:---|:---|
| **`chb01.m`** | Analysis of `chb01_03.set` data: tensor extraction, CPD, rank analysis, uniqueness verification, and visualization. |
| **`chb02.m`** | Analysis of `chb02_16+.set` data following the same pipeline. |
| **`chb03.m`** | Analysis of `chb03_01.set` data following the same pipeline (note: rank = 3 for the 10-second tensor). |

## Input Data
- Download the required `.set` files from the [CHB-MIT Database](https://physionet.org/content/chbmit/1.0.0/).
- Files needed:
  - `chb01_03.set`
  - `chb02_16+.set`
  - `chb03_01.set`

## Requirements
- MATLAB
- EEGLAB toolbox (for loading `.set` files)
- Tensor Toolbox (for tensor decomposition)

## Notes
- Tensor decomposition was repeated multiple times with different initializations to verify uniqueness.
- The low rank was determined based on core consistency and relative error plots.

