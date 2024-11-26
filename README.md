[![DOI](https://zenodo.org/badge/827257897.svg)](https://doi.org/10.5281/zenodo.12722805)

# Can we Constrain Geographical Variability in the Biological Carbon Pump's Transfer Efficiency from Observations?

![Comparison of published BCP mesopelagic transfer efficiency metrics across six
ocean sites associated with time-series programs](https://github.com/user-attachments/assets/a97488b3-3234-4993-b64c-3a38e6624162)

This repository contains the MATLAB scripts and datasets used for the data analysis and figure generation in the study: 

**"Can We Constrain Geographical Variability in the Biological Carbon Pump's Transfer Efficiency from Observations?"**

A. Rufas<sup>1</sup>, S. Khatiwala<sup>1</sup>, K. M. Bisson<sup>2,</sup><sup>3</sup>, A. P. Martin<sup>4</sup>, H. A. Bouman<sup>1</sup>

<sup>1</sup>Department of Earth Sciences, University of Oxford, Oxford, UK
<sup>2</sup>Ocean Biology and Biogeochemistry Program, NASA Headquarters, Earth Science Division, Washington, D.C., USA
<sup>3</sup>Department of Botany and Plant Pathology, Oregon State University, Corvallis, OR, USA
<sup>4</sup>National Oceanography Centre, Southampton, UK

Find the pre-print in the [Earth and Space Science Open Archive](https://essopenarchive.org/users/806280/articles/1197117-can-we-constrain-geographical-variability-in-the-biological-carbon-pump-s-transfer-efficiency-from-observations).

## Requirements

To use the content of this repository, ensure you have the following.
- [MATLAB](https://mathworks.com/products/matlab.html) version R2021a or later installed. 
- Third-party functions downloaded from [MATLAB's File Exchange](https://mathworks.com/matlabcentral/fileexchange/): `worstcase`, `MCErrorPropagation`, `swtest.m`, `FMINSEARCHBND`, `m_map`, `brewermap`, `longhurst_v4_2010`, `subaxis` and `plotBarStackGroups`. Once downloaded, please place the functions in the `./resources/external/` directory.
- MATLAB toolboxes: the [Optimization Toolbox](https://mathworks.com/products/optimization.html), necessary to run `worstcase`, the [Statistics and Machine Learning Toolbox](https://mathworks.com/products/statistics.html), necessary to run `MCErrorPropagation` and `swtest.m`, and the [Curve Fitting Toolbox](https://www.mathworks.com/products/curvefitting.html), necessary to run `fittype`.

## Repository Structure

 - `./code/`: contains the MATLAB scripts for analysing and visualising data (*provided, see "MATLAB Scripts" section*).
 - `./data/`
    - `./raw/`: raw data compiled for this study (*partly provided, see "Obtaining Raw Data" section*).
    - `./interim/`: input processed data for MATLAB scripts (*provided, see "Obtaining Interim Data" section*).
    - `./processed/`: intermediate or input data used by the MATLAB scripts for further processing (*provided*).
- `./resources/`
    - `./external/`: third-party resources for plotting and functions (*see "Requirements" section*).
    - `./internal/`: custom MATLAB functions generated specifically for calculating (*provided*).
- `./figures/`: figures generated from processed data (*provided*).

## Obtaining Raw Data

- **Particulate organic carbon (POC) flux measurements** from sediment traps and radionuclides, a compilation that we made for this study. The file, `dataset_s0_trap_and_radionuclide_compilation.xlsx`, is not available within this repository as it contains data owned by other authors that are not in a preservation repository. References for constructing this dataset are provided in the Supporting Information of our paper (Tables S1-S6), and a processed version of the dataset is available in the folder `./data/processed/pocflux_compilation.mat`.
- **Particle concentration measurements** from the Underwater Vision Profiler 5 (UVP5) downloaded from the EcoPart repository hosted by [IFREMER](https://ecopart.obs-vlfr.fr). This data are stored in the subfolder `./data/raw/UVP5/` and also in a processed format in `./data/processed/UVP5/pocflux_bisson_45sc.mat`.

## Obtaining Interim Data

- World Ocean Atlas 2023 annual climatology for **temperature**. This dataset, which I downloaded from the [NCEI NOAA website](https://www.ncei.noaa.gov/products/world-ocean-atlas), is stored in the file `temp_annual_woa23.mat`. It is required for running the Marsay et al. ([2015](https://doi.org/10.1073/pnas.1415311112)) algorithm (refer to Text S5 in the Supporting Information).
- The VHRR Pathfinder v.5.0 global 4 km monthly climatology (1985–2001) for **sea surface temperature** (SST). This dataset, which I downloaded from the [NCEI NOAA website](https://www.ncei.noaa.gov/products/avhrr-pathfinder-sst), is stored in the file `sst_pathfinder_v5.mat`. It is required for running the Henson et al. ([2012](https://doi.org/10.1029/2011GB004099)) algorithm (refer to Text S5 in the Supporting Information).
- Carr 2002 **net primary production** (NPP) monthly climatology. I constructed this climatology using the Carr [2002](https://doi.org/10.1016/S0967-0645(01)00094-7) algorithm for NPP, which incorporated inputs from (i) the SeaWiFS 9 km monthly climatology (1997-2010) for **chlorophyll concentration** and **photosynthetically available radiation** (PAR<sub>0</sub>), available for download from the [NASA website](https://oceancolor.gsfc.nasa.gov/l3/), and (ii) the AVHRR Pathfinder v.5.0 global 4 km monthly climatology (1985–2001) for SST. The resulting file, `npp_carr2002_seawifs_pathfinder.mat`, is essential for running the Henson et al. ([2012](https://doi.org/10.1029/2011GB004099)) algorithm (refer to Text S5 in the Supporting Information).
- **Euphotic layer depth** (*z<sub>eu</sub>*) monthly climatology. I constructed this climatology using two data products from the E.U. Copernicus Marine Service (CMEMS): the **diffuse attenuation coefficient at 490 nm** ([k<sub>d</sub>(490)](https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/description)) from the Copernicus-GlobColour merged product based on satellite observations, and the **mixed layer depth** ([MLD](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/services)) from the physics reanalysis product. I have used the 0.1% light depth for the calculation, as recommended by Buesseler et al. ([2020](https://www.pnas.org/doi/10.1073/pnas.1918114117)). The file, `zeu_calculated_kdcmems_mldcmems_pointonepercentpar0.mat`, is necessary to run the functions for fitting Martin's *b* and the remineralisation length scale coefficient (*z**).

The scripts that I used to read the `.nc` files downloaded from the web repositories and generate the `.mat` files presented above are available in this related [repository](https://github.com/annarufas/ocean-data-lab).

## MATLAB Scripts

The following scripts have been run in this order to analyse the data and reproduce the figures in our paper.

| Num| Script name                                  | Script action                                                     |
|----|----------------------------------------------|--------------------------------------------------------------------
| 1  | plotGlobalMapWithTimeseriesStations.m        | Creates **Figure 1**                                              |
| 2  | processPocFluxFromTrapAndRadCompilation.m    | Processes **Dataset S0**                                          |
| 3  | plotPocFluxFromTrapAndRadCompilation.m       | Creates **Figure 2**, **S1** and **S2**                           |
| 4  | processPocFluxFits.m                         | Calculates *b* for various scenarios using **Dataset S0** and generates **Dataset S1**     |
| 5  | plotPocFluxFits.m                            | Creates **Figure S3**                                                 |   
| 6  | processPocFluxFromUvp.m                      | Processes the **UVP5 dataset** downloaded from Ecopart            |
| 7  | plotPocFluxFromUvp.m                         | Creates **Figure S5**                                             | 
| 8  | findAndPlotUvpVsCompilationPocFluxMatchups.m | Creates **Figure 3**                                              |
| 9  | processBcpMetrics.m                          | Calculates *b*, *z** and T<sub>eff</sub> using the best method determined after script 7 and generates **Dataset S2** |
| 10 | plotBcpMetrics.m                             | Creates **Figure 4** and  **Figure S6** |
| 11 | examineEffectSampleSize.m                    | Analyses Dataset S0 to examine effect of sample size in the calculations of the relative error of *b* and *z** |

The following are helper functions used by the scripts above:

| Num| Script name                                    | Script action                                                     |
|----|------------------------------------------------|--------------------------------------------------------------------
| 12 | calculateBcpMetricsFromTrapAndRadCompilation.m | Steps to calculate and propagate error for *b*, *z** and T<sub>eff</sub> from Dataset S0; called by scripts 4 and 9 |
| 13 | calculateBcpMetricsFromUvp.m                   | Steps to calculate and propagate error for *b*, *z** and T<sub>eff</sub> from the UVP5 dataset; called by script 9 |
| 14 | calculateBcpMetricsFromHenson2012.m            | Steps to calculate and propagate error for *b*, *z** and T<sub>eff</sub> using the algorithm of Henson et al. ([2012](https://doi.org/10.1029/2011GB004099)); called by script 9 |
| 15 | calculateBcpMetricsFromMarsay2015.m            | Steps to calculate and propagate error for *b*, *z** and T<sub>eff</sub> using the algorithm of Marsay et al. ([2015](https://doi.org/10.1073/pnas.1415311112)); called by script 9 |
| 16 | propagateErrorWithMCforMartinbAndZstar.m       | Monte Carlo error propagation algorithm for *b* and *z** fits; called by scripts 12 and 13  |
| 17 | propagateErrorWithMCforTeff.m                  | Monte Carlo error propagation algorithm for T<sub>eff</sub>; called by scripts 12 and 13 |
| 18 | solveMartinbAndZstar.m                         | Algorithm to solve *b* and *z**; called by scripts 11 and 16        |
| 19 | samplePocFluxWithLH.m                          | Latin Hypercube sampling algorithm for POC flux profiles; called by scripts 11 and 16 |
| 20 | binPocFluxDataAtRegularDepthIntervals.m        | Low-level function, bins POC flux data by depth intervals; called by scripts 11 and 12          |
| 21 | extractDataFromZrefToZmeso.m                   | Low-level function, extracts POC flux data from the reference depth to the base of the mesopelagic zone; called by scripts 5, 11, 12 and 13 |
| 22 | extractDataFromZrefToEnd.m                     | Low-level function, extracts all POC flux data below the reference depth; called by script 5 |
| 23 | constructFilenameFitMetrics.m                  | Low-level function, constructs file names for *b* fits for different scenarios; called by scripts 5, 6 and 16 |

## Reproducibility

Our scripts showcase the application of Monte Carlo-based techniques for error propagation across diverse datasets of biological carbon pump (BCP) mesopelagic transfer efficiency metrics. Specifically, the scripts are hard-wired to handle data from our six designated study sites:
- the Hawaii Ocean Time-series (HOT) station ALOHA (HOT/ALOHA), in the subtropical NE Pacific (22.45ºN, 158ºW);
- the Bermuda Atlantic Time-Series/Oceanic Flux Program joint site (BATS/OFP), in the subtropical NW Atlantic (31.6ºN, 64.2ºW);
- the US JGOFS Equatorial Pacific process study experimental site (EqPac), in the central equatorial Pacific upwelling system (–2 to 2ºN, 140ºW);
- the Porcupine Abyssal Plain time-Series Observatory (PAP-SO), in the subpolar NE Atlantic (49.0ºN, 16.5ºW);
- Ocean Station Papa (OSP), in the HNLC region of the subpolar NE Pacific (50ºN, 145ºW), and
- the Long-Term Ecological Research (LTER) observatory HAUSGARTEN, in the polar Atlantic-Arctic boundary (79ºN, 4.0ºE).

## Acknowledgments

This work was completed as part of my PhD project at the University of Oxford under the NERC large grant COMICS (Controls over Ocean Mesopelagic Interior Carbon Storage, NE/M020835/2). I also acknowledge funding from the University of Oxford's Covid-19 Scholarship Extension Fund and Oxford's Wolfson College Covid-19 Hardship Fund.
