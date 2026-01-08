# Ice Thickness Model considering Sliding Law (ITMSL)

## Model Overview
This model is built upon the combined framework of laminar flow theory and glacier basal sliding laws to invert glacier thickness. Developed in MATLAB, it takes glacier surface DEM, surface velocity, and glacier inventory data as inputs to compute the internal thickness distribution of glaciers.

### Software Environment
- MATLAB R2018a or higher
- Required MATLAB Toolboxes:
  - Topotoolbox-master

## File Structure
```
├── velocity_to_thickness_main.m          # Main execution script
├── function/                           # Utility functions directory
│   ├── cal_basal_slid.m                 # Calculate the basal sliding
│   ├── cal_effect_k.m                   # Calculate the effective pressure coefficient
│   ├── cal_shape_factor.m               # Calculate the valley shape factor
│   ├── cal_sliding_parameter_As.m        # Calculate the sliding parameter
│   ├── dis2poly.m                      # Calculate Euclidean geometric distance
│   ├── Gantayat14_initial_icethickness.m  # Calculation of the initial thickness of 
the glacier based on Gantayat et al., 2014
│   ├── make_glacier_mask.m             # Creating Glacier Mask Files
│   └── moving_average.m                # Perform a moving average on the results
├── data/                               # Input Data
│   ├── dem.tif              
│   ├── velocity.tif
│   ├── outline.shp
│   ├── outline.dbf
│   ├── outline.prj
│   ├── outline.sbn
│   ├── outline.sbx
│   └── outline.shx
├── Previous results/             
│   ├── composite.tif           # Farinotti et al., 2019
│   ├── glabtop2.tif            # Farinotti et al., 2019
│   ├── hf.tif                 # Farinotti et al., 2019
│   ├── millan.tif             # Millan et al., 2022
│   └── oggm.tif              # Farinotti et al., 2019
├── result/                  # Utility functions directory
│   ├── itmsl_parameters.mat
│   ├── Gantayat14_ice_thickness.tif
│   ├── itmsl_basal_slid.tif
│   ├── itmsl_ice_deformation.tif
│   └── itmsl_ice_thickness.tif
└── README.md              # This file
```

## Installation Steps

1. Download
   git clone https://github.com/pxgxhcxhs/ITMSL

2. Set MATLAB Path
   - Method 1: Run in MATLAB
     ```matlab
     addpath(genpath('Your Local Path'));
     savepath; % Optional, saves path permanently
     ```
   - Method 2: Via MATLAB Menu
     - Click "Set Path" → "Add with Subfolders" → Select repository directory

## Input Data Preparation

### Required Data Files
1. Glacier Surface DEM (.tif)
   - Format: GeoTIFF
   - Requirements: Contains elevation data in meters

2. Glacier Surface Velocity (.tif)
   - Format: GeoTIFF
   - Requirements: Contains velocity data in m/year

3. Glacier Inventory File (.shp)
   - Format: Shapefile (including .shp, .shx, .dbf, etc.)
   - Requirements: Defines glacier boundary extent

### Data Consistency Requirements
- All raster data must have the same:
  - Spatial resolution (e.g., 30 m)
  - Coordinate reference system (UTM recommended)
  - Spatial extent (clipping to consistent extent recommended)

### Example Data
The `data/` directory contains test data:
- `dem.tif`: Example DEM
- `velocity.tif`: Example velocity
- `outline.shp`: Example glacier boundary

### Quick Start
Run with example data:
```matlab
% Run complete workflow
run_example('velocity_to_thickness_main.m')
```

## License

This project is licensed under the BSD-3-Clause license.

## Changelog

## Contact Information

- Author: Xiaoguang Pang
- Email: pangxg@haut.edu.cn
- Project Homepage: https://github.com/pxgxhcxhs/ITMSL
