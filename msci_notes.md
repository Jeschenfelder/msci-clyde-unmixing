# MSci 2021 - Quantitative Provenance: Chemical Fluxes across landscapes

Project working on quantifying and modelling geochemistry across river systems.  

- **Forward problem:** Using chemistry in geochemical maps to predict downstream river sediment chemistry
- **Inverse problem:** Using chemistry from river sediment samples to predict first-order geochemical maps

## Methods

### Preprocessing

#### DEM

- Using SRTMS 1 arcsecon DEM (source: https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1-arc?qt-science_center_objects=0#qt-science_center_objects)
- Downsampled to 100x100m grid (gmt grdsample)
- Projected to cylindrical equal area projection (gmt grdproject)
- Sink filling to allow correct flow-routing

#### Data and Flow

##### Setting up topography and localities

- Given "Compiled\_Clyde\_Surface" data unique sample numbers
	+ 1000 to 1010
- Filter surface data to not include estuary samples
	+ Only taking samples with lon > -4.42 (from map)
	+ Using "preprocess\_topo\_loc" script
- Localities are snapped to flow model grid in python
	+ manually checked and some nudged to fit correctly

##### Flow routing and Drainage Areas

- D8 flow routing algorithm to estimate streams
	+ Need to tweak min catchment area for ideal settings
- Snap sample data to nearest stream
	+ Need to manually check each data point to see whether correct snapping happened
	+ Use matplotlib widgets in notebook and gmt maps to correlate
	
##### Geochemical data:

- Combined Surface and CUSP data sets in one sheet
- Convert all oxides to real element concentrations
	- Constants https://www.jcu.edu.au/advanced-analytical-centre/resources/element-to-stoichiometric-oxide-conversion-factors
- Following locations needed nudging to snap:
	+ 700012: fit on main channel

##### G-BASE input

- Data from: -----
- filter through to remove bad data and pick correct sampling method for each locality using script by Alex
- Interpolate using `interpolate_gbase.gmt` script
	+ save unprojected and projected version to plot and use
- changing dampening factor 
	
### Model

- Base code for inverse from Lipp et al. 2021
	+ Modified to work for new data set and region
	+ Run on HPC to find ideal lambda for each element
- Forward model using landlab find\_drainage\_area\_and\_discharge
	+ Needs interpolated G-BASE data for region
	+ interpolation using gmt surface to 100x100m grid
	
### Forward Model Workflow

1. Interpolate G-BASE data for element of choice using `interpolate_GBASE.gmt`, move grids to `DATA/INTERPOLATED_GBASE`; MUST CHANGE OUTPUT FILENAME AND $ IN AWK STATEMENT!
2. Change element name in `forward_model.py` and run, output is saved to `DATA/FORWARDMODEL_RESULTS/`
3. Plot full results output with `diagnostics_dashboard.gmt`
	
## Results

### Forward Model

- In Mg, see most misfit in Glasgow area, draining NE of Glasgow
	+ dataset differences?
	+ lack of G-BASE data in Glasgow? misfits line up with river segment only fed by glasgow city area


## Further ideas:

- Quantify effect of seawater cation exchanges in Firth
	+ Using forward model
	+ read [Major Element Composition of Sediments in Terms of Weathering and Provenance: Implications for Crustal Recycling](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GC008758)
- Quantify effect of Glasgow on river chemistry
	+ Using forward model
	+ running with/without Glasgow input
- Improve model using machine learning
	+ read [A machine learning approach to geochemical mapping](https://www.sciencedirect.com/science/article/pii/S037567421630098X)
- Parallelising code using OpenMP/Cython


## Literature

### Cation Exchange

- Sediment-seawater cation exchange in Himalaya, Lupker et al 2016: https://esurf.copernicus.org/articles/4/675/2016/esurf-4-675-2016.pdf
	+ Measured cation exchange rates in Ganga and Brahmaputra basins
	+ Using cation exchange capacity (CEC), number of cations bound to mineral surface charges that can be reversibly exchanged
	+ CEC strongly variable in water column -> sediment sorting with depth
	+ CEC of all cations very close to Ca^2+ and Mg^2+ values
	+ dissolved Ca^2+ in ocean higher due to exchanges with Na^+ (6% flux), Mg^2+ & K^+ (5% flux)
	+ readsorbtion of Na^+ , Mg^2+ and K^+ also existent
	+ cation exchange does not signficantly impact effect of silicate weathering on carbon sequestration
	+ coarse, low surface area lead to overall low CEC

### Anthropogenic impacts

- Government guidlines on heavy metal concentrations in river sediments: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/291646/scho1108bozd-e-e.pdf




## Meeting Notes

### 8th October - Roberts

- get forward model to run first
	+ for Magnesium first
	+ recreate Alex's figures about it
	+ get code from github and 
	+ try some other elements as well

### 10th October - Lipp

- need to write Powell with iteratively printing
	+ test first with 10x10
- use find\_drainage\_area\_and\_discharge to run forward model
 

### 12th October - Lipp

- extract data misfit from minimum working example (prediction <-> observations)
- check versus G-BASE data
- port model to .py scripts
	+ generic loading in
	+ for loop through different lambda
- get G-BASE interpolation
- Try Chromium
- for HPC
	+ be on Imperial VPN
	+ username needs clearing
	+ uses pbs for jobs

### 18th October - Roberts & Lipp

- Overall forward model results promising
- interpolated G-BASE data looks spotty
	+ Try smoother interpolation by decreasing tension factor in grd surface
	+ check G-BASE points plotted to see expected interpolation
- tighter colourscale could help show variability better
	+ homogeneity can be from lack of changes too
- RMS should be around 0.1 if model predicts well
	+ R^2 of less use as more variabilty is to be expected
- In cross plots, plot CUSP and surface data separately
- [ ] Create script for big diagnostics plot -> can be generated for each element; make two scripts? (GMT 1-3; Python 4-6)
	1. 	raw distribution
	2. 	interpolated G-BASE data
	3. 	results with observations above
	4. 	cross plot; showing rms and misfit histogram
	5. 	map of misfits
	6. 	river profile
- Should get Glasgow G-BASE data soon
	+ not pristine (used from sludge)
	+ most likely polluted
- Could use geology to find natural baseline
	+ average concentrations from underlying lithology before interpolation 
- [ ] Overlying geologic boundaries on interpolation to see larger structures
- Next meeting: Thursday for next week plan
- gmt tip, run basemap last to avoid artefacts on sides (need to check if this works in modern version)

### 19th October - Lipp

- Run on HPC
	+ add array_index
		* use counter to have each node run just one inversion
	+ create job scripts
		* starting at 2gb memory
		* use 2 CPU to avoid crashed
	+ also output model roughness and data misfit for each lambda
- Setup HPC
	+ using anaconda to install landlab
	+ regular 
- Submit jobs with `qsub`
	+ check jobs with `qstat`
	+ array jobs show up wiht [] after ID
	+ remove faulty jobs with `qdel`
- Try higher resolutions
	+ keep squares
	+ feel free to go high with res
	+ try ~60x50 first (keep square)
- Get landlab working in anaconda
	+ set up anaconda (help in presentation)
- using gmt classic more helpful for dashboards
	+ sent a snippet to get started
	+  -O (overwrite) -K (don't finalize) flags important
	+  

## Weekly Plan/Progress

### 20th - 26th September:

- Meeting to talk about first project plans
- Reading papers
	+ Lipp 2020
	+ Lipp 2021
- Check out published code for models
- Start learning gmt again

### 27th - 3rd October:

- Create intro plots of region and write intro presentation
- Travel to UK 28th October all day
- Intro presentation at EGG-fest (Geochemistry and Geodynamics research group meeting) 29th October
- Getting full code for inverse model

### 4th - 10th October:

- Start adapting inverse code to new region
- Fitting localities to drainage network
- Intro sessions/Welcome Fair

### 11th - 17th October:

- Finishing setting up inverse code in jupyter notebook
	+ Getting small example (Mg, 10x10 grid, lambda = 2) to run
- Writing forward model code from inverse
	+ preprocessing same
	+ running forward model with `landlab find_drainage_area_and_discharge()`
	+ plotting results, misfit by locality and cross plot
- Start converting inverse code from jupyter notebook to script
	+ Write to have variable lambda
	+ Parallelisable for use on HPC
	+ Output each raster as `.asc` file

### 18th - 24th October:

- Finished inverse script for HPC
- Start on diagnostics plotting scripts
	+ Decided to do maps in gmt and rest in 