# MSci 2021 - Quantitative Provenance: Chemical Fluxes across landscapes

Project working on quantifying and modelling geochemistry across river systems.  

- **Forward problem:** Using chemistry in geochemical maps to predict downstream river sediment chemistry
- **Inverse problem:** Using chemistry from river sediment samples to predict first-order geochemical maps

## Methods

### Preprocessing

#### DEM

- Using SRTMS 1 arcsecond DEM (source: https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1-arc?qt-science_center_objects=0#qt-science_center_objects)
- Downsampled to 100x100m grid (gmt grdsample)
- Projected to cylindrical equal area projection (gmt grdproject)
- Sink filling to allow correct flow-routing

#### Data and Flow

#### Setting up topography and localities

- Given "Compiled\_Clyde\_Surface" data unique sample numbers
	+ 1000 to 1010
- Filter surface data to not include estuary samples
	+ Only taking samples with lon > -4.42 (from map)
	+ Using "preprocess\_topo\_loc" script
- Localities are snapped to flow model grid in python
	+ manually checked and some nudged to fit correctly

#### Flow routing and Drainage Areas

- D8 flow routing algorithm to estimate streams
	+ Need to tweak min catchment area for ideal settings
- Snap sample data to nearest stream
	+ Need to manually check each data point to see whether correct snapping happened
	+ Use matplotlib widgets in notebook and gmt maps to correlate
	
#### Geochemical data:

- Combined Upper Clyde River Sediment Survey and Estuarine Geochemical (surface) data sets in one sheet
- Choosing correct concentrations for Upper Clyde Survey (for duplicates): from Paul Everett (BGS)
	+ CaO, TiO2, Fe2O3 and Cl from WDT
	+ Ba from ED
	+ if Fe: 20+%, Ti: 10+%, Ca: 30+%, Cl: 1%, Ba: 0.5+% --> use WDM data 
- Convert all oxides to real element concentrations
	- Constants https://www.jcu.edu.au/advanced-analytical-centre/resources/element-to-stoichiometric-oxide-conversion-factors
- Following locations needed nudging to snap:
	+ 700012: fit on main channel

#### G-BASE input

- Data from:G-BASE geochemical survey 
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
	- `1.sample_no 2.easting 3.northing 4.Ag\_out 5.Ag\_meth 6.Al\_out 7.Al\_meth 8.As\_out 9.As\_meth 10.B\_out 11.B\_meth 12.Ba\_out` 
`13.Ba\_meth 14.Be\_out 15.Be\_meth 16.Bi\_out 17.Bi\_meth 18.Br\_out 19.Br\_meth 20.Ca\_out 21.Ca\_meth 22.Cd\_out 23.Cd\_meth 24.Ce\_out` 
`25.Ce\_meth 26.Cl\_out 27.Cl\_meth 28.Co\_out 29.Co\_meth 30.Cr\_out 31.Cr\_meth 32.Cs\_out 33.Cs\_meth 34.Cu\_out 35.Cu\_meth 36.Fe\_out`
`37.Fe\_meth 38.Ga\_out 39.Ga\_meth 40.Ge\_out 41.Ge\_meth 42.Hf\_out 43.Hf\_meth 44.Hg\_out 45.Hg\_meth 46.I\_out 47.I\_meth 48.In\_out`
`49.In\_meth 50.K\_out 51.K\_meth 52.La\_out 53.La\_meth 54.Li\_out 55.Li\_meth 56.Mg\_out 57.Mg\_meth 58.Mn\_out`
`59.Mn\_meth 60.Mo\_out 61.Mo\_meth 62.Na\_out 63.Na\_meth 64.Nb\_out 65.Nb\_meth 66.Nd\_out 67.Nd\_meth 68.Ni\_out 69.Ni\_meth 70.P\_out 71.P\_meth 72.Pb\_out`
`73.Pb\_meth 74.Rb\_out 75.Rb\_meth 76.S\_out 77.S\_meth 78.Sb\_out 79.Sb\_meth 80.Sc\_out 81.Sc\_meth 82.Se\_out 83.Se\_meth 84.Si\_out`
`85.Si\_meth 86.Sm\_out 87.Sm\_meth 88.Sn\_out 89.Sn\_meth 90.Sr\_out 91.Sr\_meth 92.Ta\_out 93.Ta\_meth 94.Te\_out 95.Te\_meth 96.Th\_out`
`97.Th\_meth 98.Ti\_out 99.Ti\_meth 100.Tl\_out 101.Tl\_meth 102.U\_out 103.U\_meth 104.V\_out 105.V\_meth 106.W\_out 107.W\_meth 108.Y\_out` 
`109.Y\_meth 110.Yb\_out 111.Yb\_meth 112.Zn\_out 113.Zn\_meth 114.Zr\_out 115.Zr\_meth`
2. Change element name in `forward_model.py` and run, output is saved to `DATA/FORWARDMODEL_RESULTS/`
3. Get R2 and RMS misfit stats from `forward_stats.py`, output is saved to `DATA/FORWARDMODEL_RESULTS/` 
4. Plot full results output with `diagnostics_dashboard.gmt` need to change $ in AWK STATEMENTS! ADD COEFFICIENT AND CHANGE TO PPM IF WORKING ON OXIDE/MAJOR ELEMENT!
	- Upper River: ```1.SAMPLE\_No	2.Easting	3.Northing	4.SAMPLE_TYPE	5.METHOD   6.Na2O	7.MgO	8.Al2O3	9.SiO2	10.P2O5	11.SO3	12.CaO```	  
	```13.TiO2	14.Fe2O3 15.Ba	16.Cl	17.K2O	18.CaO	19.TiO2	20.MnO 21.Fe2O3	22.S	 23.Cl	24.Sc```
	```25.V 26.Cr	27.Co	28.Ni	29.Cu	30.Zn	31.Ga	32.Ge 33.As	34.Se	35.Br	36.Rb```
	```37.Sr	38.Y 39.Zr	40.Nb	41.Mo	42.Nd	43.Sm	44.Yb 45.Hf	46.Ta	47.W	 48.Tl```
	```49.Pb	50.Bi	51.Th	52.U	 53.Ag	54.Cd	55.In	56.Sn 57.Sb	58.Te	59.I 60.Cs	61.Ba	62.La	63.Ce```
	- Clyde Estuary: ```1.SAMPLE_No	2.SAMPLE_NAME	3.CRUISE	4.Lab_Number	5.Sample	6.Site	7.Easting	8.Northing	9.LOI	10.Hg	11.TOC	12.Na2O```
	```13.MgO	14.Al2O3	15.SiO2	16.P2O5	17.K2O	18.CaO	19.TiO2	20.MnO 21.Fe2O3	22.S	23.Cl	24.Sc```
	```25.V	26.Cr	27.Co	28.Ni	29.Cu	30.Zn	31.Ga	32.Ge 33.As	34.Se	35.Br	36.Rb```
	```37.Sr	38.Y	 39.Zr	40.Nb	41.Mo	42.Ag	43.Cd	44.In 45.Sn	46.Sb	47.Te	48.I```
	```49.Cs	50.Ba	51.La	52.Ce	53.Nd	54.Sm	55.Yb	56.Hf 57.Ta	58.W	59.Tl	60.Pb	61.Bi	62.Th	63.U```
	
### Inverse Model Workflow

1. Check sample data if element has enough datapoints
2. Run inverse model on HPC for smoothing -2 to 2 (41 steps); `inverse_model_HPC.py`
	- Copy job script from template and change to correct element
	- Outputs: `$elem_$lambda_roughness_misfit.txt`; `$elem_$lambda_inverse_output.asc.npy`
3. Once all done, `scp` all `$elem_*_roughness_misfit.txt` and plot tradeoff curve using `plotting_tradeoff_curve.py`
	- Pick best inverse at bend of elbow
4. `scp` `$elem_$lambda_inverse_output.asc.npy` for best inverse (from 3); turn into readable file using `inverse_npy_to_asc.py`
5. Plot inverse dashboard using ??? 
	- Remember to create G-BASE masked for element using `mask_gbase.gmt`
	- Use `evaluate_inverse.py` to create block averaged G-BASE

	
## Results

### Forward Model

- In Mg, see most misfit in Glasgow area, draining NE of Glasgow
	+ dataset differences?
	+ lack of G-BASE data in Glasgow? misfits line up with river segment only fed by glasgow city area
- River profile for Mg very uniform
- Group of estuarine samples always misfitting
	+ Only points in that dataset not taken from a boat
	+ Different methodology?
	+ Try excluding them from inverse?

### Inverse Model

- 

## Project Plan

### Figures

### Outline

- Introduction
- Methods
- Results
- Discussion
- Conclusion

## Further ideas:

- ~~Quantify effect of seawater cation exchanges in Firth~~
	+ ~~Using forward model~~
	+ ~~read [Major Element Composition of Sediments in Terms of Weathering and Provenance: Implications for Crustal Recycling](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019GC008758)~~
- Quantify effect of Glasgow on river chemistry
	+ Using forward model
	+ running with/without Glasgow input
- Improve model using machine learning
	+ read [A machine learning approach to geochemical mapping](https://www.sciencedirect.com/science/article/pii/S037567421630098X)
- Parallelising code using OpenMP/Cython
- Adding point sources!
- Correlate predicted contaminant sources to industrial areas in basin

## Literature

### Cation Exchange - NOT PURSUED

- Sediment-seawater cation exchange in Himalaya, Lupker et al 2016: https://esurf.copernicus.org/articles/4/675/2016/esurf-4-675-2016.pdf
	+ Measured cation exchange rates in Ganga and Brahmaputra basins
	+ Using cation exchange capacity (CEC), number of cations bound to mineral surface charges that can be reversibly exchanged
	+ CEC strongly variable in water column -> sediment sorting with depth
	+ CEC of all cations very close to Ca^2+ and Mg^2+ values
	+ dissolved Ca^2+ in ocean higher due to exchanges with Na^+ (6% flux), Mg^2+ & K^+ (5% flux)
	+ readsorbtion of Na^+ , Mg^2+ and K^+ also existent
	+ cation exchange does not signficantly impact effect of silicate weathering on carbon sequestration
	+ coarse, low surface area lead to overall low CEC

### Anthropogenic impacts - MAIN IDEA

- Government guidlines on heavy metal concentrations in river sediments: https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/291646/scho1108bozd-e-e.pdf
	+ Found increased levels of Cd, Pb, Cu, Zn, As in rivers across England and Wales
	+ Discharge during active periods higher, but reservoirs still exist in river beds
	+ Hazards: mainly due to Re-suspension and dispersion during floods
		* Harm to aquatic life
		* Contaminate soil and agriculture around
	+ Long history of mining in England and Wales -> often using rivers to sort ore
	+ Key long=term impact is contamination of fine-grained suspended sediments
		* especially Cd, Cu, Pb, Zn
		* i.e. Devon and Cornwall show significant pollution
	+ Main contamination associated with late 19th century mining (via deposition depth)
	+ Contamination in some cases leads to unique ecosystems forming (i.e. Pennine Orefield)
	+ Flodding can distribute vast amounts of heavy metals over an area
		* i.e. Autum 2000 flooding causing dispersial of high Pb soils up to 80km away (1000ppm found)
	+ Climate change related increased frequency in flooding making this a bigger concern
	+ Rivers as temporary storage of contaminant, floodplains, wtlands etc as long-term sinks
		* 10s-100s kyrs in floodplains
		* storage in rivers can exceed 100 yrs (especially in estuarine sediments)
	+ Safe levels of pollution in sediments are not yet set
	+ Debate about toxicity in soils
	+ TEL and PEL as interim solutions based on Canadian environmental standards
		* TEL: threshold effect level; concentration below no significant hazard to aquatic organisms is posed
		* PEL: predicted effect level; concentrations associated with adverse biological effects
		* see Table 2.4 for proposed levels (page 15)
	+ sediment borne contaminants affect flora, fauna and humans
	+ Studies show strong uptake of Cd, Pb and Zn in crops from soil
		* Excess metals can build up in the plants or excluded via the routs
		* Grazing stock often ingest that plant material or sediment -> risk of exposure through grazing
	+ Risks to human health possible:
		* 

### Point sources

- Via concentration:
- Via mass: 



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
- [x] Create script for big diagnostics plot -> can be generated for each element; make two scripts? (GMT 1-3; Python 4-6)
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
- [x] Overlying geologic boundaries on interpolation to see larger structures
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
	
### 26th October - Lipp

- Suspect NaN result bug comes from `data_misfit` function
- Once inverse runs on HPC should focus on getting diagnostics for forward up
- Low variability in Mg likely due to mainly sedimentary rocks in region
	+ Can try forward model on elements with highest variability in observations next
- Given me the scripts to create synthetic data for inverse, should use eggbox first

### 2nd November - Roberts; Lipp

- Should run synthetics for inverse; sinusoids 5,10,15,20,25 km
- Can use PCA to collate different element outputs
- Storyline of final writeup:
	+ showing performance of mixing/unmixing model in new region
	+ predicting pollution sources from model outputs
	+ discuss natural vs human impact
	+ discuss data availability/quality
- Need to make sure to talk about data in project
	+ Issues with estuary samples, make clear which ones used (shown in all Figures)
	+ Don't know exact uncertainty, can use Lipp 2020 as approximation
- Figures for writeup:
	1. Intro to region:
		+ Topography and Drainage Network
		+ Observation localities
		+ Geological Map
		+ G-BASE sample localities
	2. 4 interesting elements:
		+ 2 good ones (i.e. Mg)
		+ 2 'bad' ones (i.e. Metals) -> Human influence
- Before running more inverse models, change to new model in C -> much faster
- Changes to forward dashboards (Appendix):
	+ letters without ')'
	+ get rid of unused observations
	+ fill profile points by misfit
	+ keep cross plot range constant at 0/5
	+ don't use red for Glasgow -> doesn't print well
	+ Add PEL/TEL levels to profile (dashed lines)
- ultimately should show inverse as enrichment relative to G-BASE -> highlight possible human influences; or to PEL and TEL
	+ can check against industrial areas
	+ baseline already includes some pollution most likely
	
### 8th November - Roberts

- Should start synthetic run
- Inverse dashboards:
	+ use G-BASE points as comparison, not interpolation!
	+ add cross plot G-BASE vs inverse
	+ map of G-BASE vs inverse misfits
	+ inverse/PEL;TEL
- Run more elements in inverse and forward

### 9th November - Lipp

-  Synthetic range might be too large, should be fine
-  inversion dashboards to set up
-  might as well run all elements
	+  pick some elements of interest first
	+  anthropogenic, high variance
	+  check with G-BASE outputs
- ask Gareth about literature for environmental geochemistry
	+ outcomes of sediment pollution
	+ check library textbooks for environmental geochemistry
- can average G-BASE onto grid resolution of inverse output
	+ sent script to work off of
	+ easier way to do is just find inverse value at G-BASE
- use ratio of G-BASE to inverse (in log?)

### 16th November - Roberts

- Run 'wrong' data through to prove that doesn't fit well
	+ Use Wales data
	+ Or could use randomised G-BASE from this area
	+ Ideally both
- Could focus on just three elements
	+ high, mid, low concentrations
	+ Split Fig 4 up into two
		* a: G-BASE, river; obs v pred; profile
		* b: cross plot of elements together
	+ Add Wales and random data as proof that it works
		* show randomized G-BASE, Wales G-BASE, 1:1 for each underneath
- Lead people down the road slowly
- Show goodness of fit first then results (inverse) for best fitting element
	+ 1:1 plot of best fitting model to observation (a)
	+ Output of inversion raw (b)
	+ comparison between predicted G-BASE and real G-BASE (c)
- Prediction of source comp for 6 elements
- Uses figure:
	+ PEL/TEL from forward
	+ PEL/TEL on inversion
- In discussion about point sources:
	+ Analyse what size of source needs what concentrations to be detectable
	+ When do you reach the point that it is not detectable anymore?
	+ Contours of different initial concentrations
	+ Upstream area against concentration
	+ Equation 2 from Lipp 2021 (Cinit*Ap_source)/Aupstream
	+ Dimensionless solution (?)
	+ Add dashed lines showing detection limits
	+ Try this in forward model if time permits


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
- Ran first inverse (Mg) suite on HPC
	+ Bug with final misfit, needs fixing
		* Found workaround on local version, needs testing on HPC version before more runs
	+ Found bug with how results got saved
		* Fixed by adding result to ``RasterModelGrid` first, then saving output as `.asc`
- Create plot for inverse results

### 25th - 31st October:

- Fixed inverse model for HPC
	+ bug was how initial NaN values were handled, using `np.nansum` and `np.nanmean` now in `initiate_inverse_smart`
- Running inverse on HPC for elements Mg, Co, Cr, K, Pb, Sr
	+ using lambda exponent values of -2 to 2
- Finishing diagnostics dashboard for forward model
- Need to figure out way to distinguish data sets in plots
	+ Main issue is they're handled the same in python and come as same output
	+ Save in separate files from get go? --> easy for misfit data
		* possible for profile?
		* Could go by distance in profile?
	+ Or doable in gmt? --> avoids cluttering up result files
- Run forward model for Pb, Sn, Sb, Zr
	+ no G-BASE data for Br, Hf --> cannot run forward
	+ need to rerun Sb as interpolation included negative concentrations
	
### 1st - 7th November:

- finished initial forward model runs
- working on inversion results -> starting dashboard
- inversion concentrations seem wrong
	+ inverse gives results in ln, need to convert to log10
	+ now mapped correctly

### 8th - 14th November:

- Run more inverse
	+ synthetics: 5;10;15;20;25km eggbox
	+ elements: Br; Sn; Zn; Zr; Cu
- Inverse predicts high copper pollution NW of Glasgow
	+ Could be due to point sources upstream Clyde not in Tributary -> adding point sources?
	+ Evidence of copper used in shipbuilding and other industries -> need to research this more
- Prepare and give short presentation on progress to research group

### 15th - 21st November:
- Create figure plan
- Create inverse dashboards

### 22nd - 28th November:
- Running forward model with randomised G-BASE and wales G-BASE data
- Create first draft of write-up figures

