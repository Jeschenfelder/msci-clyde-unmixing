# Heavy Metal Concentrations Clyde - Supplementary Material

Presented here are minimum working examples for the forward and inverse methodology. Some steps of pre- and post-processing are skipped for convenience to the reader. The aim of the materials is to give a general overview of the workflow for both models.

## File Structure:
```
supp_material/
├── Data
│   ├── converted_chem_data.csv
│   ├── filled_topography.npy
│   ├── filtered_sample_loc.dat
│   ├── gmt.history
│   ├── landmask.nc
│   ├── loc_areas.npy
│   ├── Mg_gbase_log.nc
│   └── Topography_100m.nc
├── Forward_results
│   ├── Mg_forward_result.asc
│   ├── Mg_obs_profile.txt
│   ├── Mg_obs_v_pred.txt
│   └── Mg_pred_profile.txt
├── forward_min_working_example.ipynb
├── inverse_min_working_example.ipynb
└── README.md
```

## Forward Model:

The forward model is presented in `forward_min_working_example.ipynb`. The notebook shows how to set up the model grid in LandLab and run the forward model based on an interpolated input grid. It also includes the extraction of the predicted downstream geochemistry along the main channel and basic visualisation examples.

## Inverse Model:

The inverse model is presented in `inverse_min_working_example.ipynb`. This notebook covers how to set up the model grid in Landlab, generate the active drainage network and run the inversion scheme. It is set up to run on a resolution of 10 x 10 node grid to speed it up. The inversion should take around 4 minutes to complete.


## System Requirements:
To view the materials, Python 3.8.10 or higher and jupyter 1.0.0 or higher.
To run the code the following additional packages are needed:

- landlab 2.3.0
- matplotlib 3.4.3
- netCDF4 1.5.7
- numpy 1.21.1
- pandas 1.3.2
- scipy 1.7.1


