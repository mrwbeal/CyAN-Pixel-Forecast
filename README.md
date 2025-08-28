# CyAN Pixel Forecast

Welcome to the **CyAN Pixel Forecast** project!  
This repository contains files and scripts used in the development of a pixel-scale forecasting model for Florida lakes.  

It also supports the manuscript: *Pixel-scale satellite forecasting of cyanobacteria in Florida lakes.*  

 **Note:** This repository does not contain files to run an operational forecast. It documents development and evaluation only.  

---

## Table of Contents
- [Overview](#overview)
- [File and Directory Descriptions](#file-and-directory-descriptions)
  - [Model Files](#model-files)
  - [Data Processing](#data-processing)
  - [Lead Time Evaluation](#lead-time-evaluation)
  - [Performance & Evaluation](#performance--evaluation)
  - [Other](#other)
- [Getting Started](#getting-started)
- [Author](#author)

---

## Overview
This project develops and evaluates machine learning models for pixel-scale cyanobacteria forecasting in Florida lakes.  
It includes code for data processing, model training, evaluation, and performance assessment.

---

## File and Directory Descriptions

### Model Files
- **`cyanpixelforecast_continuous_RF.R`**  
  Fold-forward evaluation of a parallel Random Forest model (`ranger`).
- **`cyanpixelforecast_continuous_ConvLSTM.R`**  
  Fold-forward evaluation of a ConvLSTM model (`keras`). Includes a data generator to handle large inputs.  
  - **`process_train_tst_convLSTM.R`**: Prepares `florida_data_spatial.csv` into tensors for ConvLSTM.
- **`cyanpixelforecast_continuous_BayesMLR.R`**  
  Fold-forward evaluation of a Bayesian MLR model (`brms`).
- **`inla_pixel_florida.R`**  
  Model code for lake-scale comparison forecast.

---

### Data Processing
- **`cyan_pixel_processing_conus.R`**  
  Scales and extracts CI data from [CyAN](https://oceancolor.gsfc.nasa.gov/about/projects/cyan/).
- **`gridMet_access.R`**  
  Accesses [GridMET](https://www.climatologylab.org/gridmet.html) THREDDS data for surface radiation and wind speed.
- **`prism_download.R`**  
  Downloads daily [PRISM](https://prism.oregonstate.edu) precipitation and temperature data.
- **`prism_pixel_processing_conus.R`**  
  Extracts pixel-scale temperature & precipitation (precip not used in model, but capability included).
- **`prism_processing_conus.R`**  
  Extracts mean temperature & precipitation values at lake scale (used in lake-scale model).
- **`prism_watershed_extract_allwatersheds.R`**  
  Extracts total weekly precipitation over lake watersheds.
- **`SWmodel_pixel_processing_conus.R`**  
  Automates surface water predictions for pixel-level extraction (adapted from Hannah Ferriby’s model).
- **`join_pixel_datasets.R`**  
  Joins predictor datasets for Florida pixel-level modeling (wind speed, surface radiance, surface water temp, precip).

---

### Lead Time Evaluation
- **`cyanpixelforecast_continuous_RF_increasedleads.R`**  
  Fold-forward Random Forest forecasts for 1, 2, and 4-week leads.
- **`Concatenate_Leads.R`**  
  Concatenates individual fold files from each lead time.

---

### Performance & Evaluation
- **`RPSS_CRPSS.R`**  
  Calculates pixel-scale RPSS and CRPSS (using WHO Alert Level 1 & 2 thresholds).
- **`ForecastPerformance.ipynb`**  
  Calculates categorical & deterministic skill scores.
- **`ForecastErrors.R`**  
  Examines errors above the 90th percentile and related predictor values.
- **`Lake_Pixel_Comparison.R`**  
  Compares performance scores between lake- and pixel-scale forecasts.

---

### Other
- **`SpatialAutocor.R`**  
  Calculates spatial correlation among pixels in each lake.
- **`TemporalAutocor.R`**  
  Calculates temporal correlation among pixels in each lake.
- **`create_pixel_shape.R`** *(optional)*  
  Creates a shapefile of Florida lake pixels.
- **`Manuscript_Figures.R`**  
  Generates figures for the manuscript.  
- **`README.md`**  
  This file.

---
### Data
- All data is available online, we have intentionally constructed this model with readily available datasets, that are updated at a weekly timestep.
- Other data needed:
  -  `MERIS_OLCI_Lakes`: the shapefile defining MERIS and OLCI resolvable lakes. 

---

## Getting Started
To reproduce parts of the workflow:
1. Clone this repository.  
2. Install [R](https://cran.r-project.org/) (≥4.0) and required modeling packages (`ranger`, `brms`, `keras`, `INLA`).  
3. (Optional) Install Python ≥3.8 with TensorFlow/Keras for ConvLSTM.  
4. Run preprocessing scripts in **Data Processing** to generate inputs.  
5. Run a model script from **Model Files**.  

**Note:** many of the data processing and model scripts require greater memory to load data into R, and greater computing power to process data in a reasonable timeframe

---

## Repository Author
- Max Beal
