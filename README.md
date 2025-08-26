# Pixel Forecast

Welcome to the **CyAN Pixel Forecast** project! This repository contains files and directories used in the development of a pixel-scale forecasting model for Florida lakes. It also supports the manuscript: Pixel-scale satellite forecasting of cyanobacteria in Florida lakes.

Author: Max Beal

Note this repository does not contain files to run an operational forecast, but rather documents development and evaluation.

## File and Directory Descriptions

### Model Files
#### `cyanpixelforecast_continuous_RF.R`
Contains model development code for a fold forward evaluation of a parallel random forest model implemented with the ranger package.

#### `cyanpixelforecast_continuous_ConvLSTM.R`
Contains model development code for a fold forward evaluation of a Convolutional LSTM model implemented with keras. Data generator code is included as well to pull tensors into memory without crashing R. 
##### `process_train_tst_convLSTM.R`
This code is needed to take the input data file "florida_data_spatial.csv" and transform into tensors.

#### `cyanpixelforecast_continuous_BayesMLR.R`
Contains model development code for a fold forward evaluation of a Bayesian MLR model implemented with the BRMS package.

#### `inla_pixel_florida.R`
Model code for the lake-scale comparison forecast.

### Increased Lead times
#### `cyanpixelforecast_continuous_RF_increasedleads.R`
Contains same development code as "cyanpixelforecast_continuous_RF.R" except develops and outputs fold foreward forecasts for 1, 2, and 4-week leads.

#### `Concatenate_Leads.R`
Concatenates the individual fold files from each lead time.

### Performance
#### `RPSS_CRPSS.R`
Takes forecast output files and calculates pixel-scale RPSS and CRPSS. RPSS is calculated using the WHO Alert Level 1 and 2 thresholds.

#### `ForecastPerformance.R`
Takes forecast output files and calculates pixel-scale categorical and determinstic skill scores.

### `README.md`
This file. It provides an overview of the project and its structure.

### Other
#### `SpatialAutocor.R`
Calcualtes spatial correlation among pixels in each lake.

#### `TemporalAutocor.R`
Calcualtes temporal correlation among pixels in each lake.
---
