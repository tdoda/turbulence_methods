# Turbulence Methods

Repository combining different methods to estimate turbulence quantities in lakes and oceans.

Link to the remote repository: https://github.com/tdoda/turbulence_methods.git

## Introduction

The aim of this project is to compare different in-situ methods for the estimation of turbulence in lakes and oceans. 

<font color='red'>*Include here a description of the project, this could include:*

- *Research objective*
- *Location*
- *Funding*
- *References*
- *Image of setup*

</font>

## Installation on a local repository 

1. To work with the project locally, download all the repository or use `git` to clone the project with the command in GitBash:

    `git clone https://github.com/tdoda/turbulence_methods.git`

    The repository will be copied to your current working directory. Note that you need to have [git](https://git-scm.com/downloads) installed in order to successfully clone the repository.

2. Use Python 3 and install the requirements with:

    `pip install -r requirements.txt`

    The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual). 

## Working with the project online without any installation (Renku)

To run **Python** scripts on the browser (without any installation), you can use the Renku
platform. Note that this approach **does not** work for Matlab scripts, that can only be ran [locally](#installation-on-a-local-repository) with Matlab installed.

Two options with Renku:
- Directly run the created session by clicking here: [![launch - renku](https://renkulab.io/renku-badge.svg)](https://renkulab.io/p/tomy.doda/turbulence-methods/sessions/01KG7ZVW6T5HWJCH73HNJFBCB9/start)
- Look for the project "Turbulence_Methods" on https://renkulab.io/ and start a new session in the `Sessions` tab. This will start an interactive environment right in your browser.

## Sensors
<font color='red'>*Include here a description of the sensors used in your project:*

- **Brand, Model & SN**
- **System integration**
- **Accuracy**
- **Setup**

</font>

## Usage
### Processing of microstructure profiles
<font color='red'>*Procedure to process microCTD and VMP profiles*</font>

## Example data

<font color='red'>*Add information the example datasets (geospatial location, temporal coverage).*</font>

## Structure of the repository 

### Folder `data`

Turbulence data organized in processing levels:
- **Level 0**: Raw data collected from the different sensors.
- **Level 1**: Raw data stored as NetCDF files where attributes (such as sensors used, units, description of data, etc.) are added to the data. 
Quality assurance is performed on the data and QA masks are generated.
- **Level 2**: Processed data, this could include calculated parameters, transformed units, resampled or gridded data.

### Folder `scripts`
Scripts used to read the data and perform further analysis, organized in the following subfolders:
- `analyze_VMP_microCTD`: contains Matlab script `analyze_profiles.m` used to process *.P datafiles from Rockland Scientific microstructure profilers (microCTD, VMP).
- `figures`: Matlab scripts used to make figures.
- `general_functions`: general Python functions called by other scripts
- `microstructure`: Matlab functions used to process microstructure profiler data.
- `odas_v4.4`: ODAS Matlab functions provided by Rockland Scientific. 

### Folder `references`
Relevant scientific papers about turbulence, organized in topic-based subfolders. 

### Folder `notes`
Anything to share (ideas, summaries, reports of analyses):
- `Turbulence_methods_manual.docx/.pdf`: manual describing the data processing for each turbulence estimation method.

 ## Future updates

<font color='red'>*Add to-do list here*</font>

## Collaborators

- **Data collection**: Tomy Doda, Haydn Herrema, Damien Bouffard, Bieito Fernández Castro,  Oscar Sepúlveda Steiner
- **Scripts and analysis**: Tomy Doda, Haydn Herrema, Oscar Sepúlveda Steiner, Bieito Fernández Castro, Sebastiano Piccolroaz, Damien Bouffard

## Contact
- [Tomy Doda](mailto:tomy.doda@unil.ch)
- [Haydn Herrema](mailto:haydn.herrema@eawag.ch)
- [Damien Bouffard](mailto:damien.bouffard@eawag.ch)


