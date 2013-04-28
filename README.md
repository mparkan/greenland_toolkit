greenland_toolkit
=================

The Greenland toolkit is a set of R scripts for exploratory analysis of in situ and remotely sensed data. 
It is primarly intented for the study of coastal areas in Greenland, but may be extended to other use cases.

It includes scripts to:
* extract and organize weather data from the NOAA National Climatic Data Center (http://www.ncdc.noaa.gov/cdo-web/)
* extract and organize precise ephemeris data from the JPL Horizons service (http://ssd.jpl.nasa.gov/?horizons)
* extract and organize summary statistics from the JPL GHRSST Level 4 AVHRR_OI Global Blended Sea Surface Temperature Analysis 
(http://podaac.jpl.nasa.gov/dataset/NCDC-L4LRblend-GLOB-AVHRR_OI)
* create monthly temperature models by location using support vector regression
* map weather stations

Software requirements
----------------------------------------------------------

* the latest release of R, available at http://www.r-project.org/
* the latest release of RStudio (IDE for R), available at http://www.rstudio.com/
* a text editor, e.g. Notepad++ (http://notepad-plus-plus.org/), SciTE (http://www.scintilla.org/SciTE.html)
* a ftp client, e.g. Filezilla (https://filezilla-project.org/), CyberDuck (http://cyberduck.ch/)
* a desktop GIS, e.g. Quantum GIS (http://www.qgis.org/), uDig (http://udig.refractions.net/), OpenJump (http://www.openjump.org/)

Workflow overview
----------------------------------------------------------

**PLEASE CHECK THE DOCUMENTATION FOR A DETAILED EXPLANATION OF THE WORKFLOW** 

Part 1 - Downloading, organizing and preparing the data sets for use with R

1. Prepare the raw data folder layout.
2. Download all the data sets and place them in their specific folders.

Part 2 - Feature construction

1. Prepare the processed data folder layout.
2. Run 'extract_ephemeris.r' to extract daily daylight duration at selected locations (detailed instructions are provided in the script).
3. Run 'extract_weather.r' to extract daily weather features at selected locations (detailed instructions are provided in the script). 
4. Run 'extract_avhrr.r' to extract daily remote sensing features at selected locations (detailed instructions are provided in the script). 
5. Run 'merge_features.r' to merge ephemeris, weather and remote sensing features (detailed instructions are provided in the script).
   
Part 3 - Data inventory and pre-analysis
   
1. Run 'locate_stations.r' to check weather data availability and compute distance matrix.
2. Run 'greenland_map.r' to plot a map of weather stations in Greenland.


Part 4 - Modelling and prediction

1. Run 'predict_self.r' to build training/test sets, create monthly SVR and linear regression models and make predictions. 
