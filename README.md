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


## Working with the project on the Renku platform (online)

The simplest way to start your project is right from the Renku
platform - just click on the `Sessions` tab and start a new session.
This will start an interactive environment right in your browser.

## Installation on a local repository 

To work with the project anywhere outside the Renku platform,
click the `Settings` tab where you will find the
git repo URLs - use `git` to clone the project on whichever machine you want.

**You need to have [git](https://git-scm.com/downloads) and [git-lfs](https://git-lfs.github.com/) installed in order to successfully clone the repository.**

- Clone the repository to your local machine using the command:

<span style="color:red">*REPLACE WITH LINK TO YOUR REPOSITORY*</span>

 `git clone 
 
 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual). 

## Sensors
<font color='red'>*Include here a description of the sensors used in your project:*

- **Brand, Model & SN**
- **System integration**
- **Accuracy**
- **Setup**

</font>

## Geospatial information

<font color='red'>*Add information about the coordinates of the profiles.*</font>

## Temporal coverage

<font color='red'>*Add information about the field campaigns.*</font>


## Structure of the repository 

### Folder `data`
Turbulence data organized in processing levels.
    - **Level 0**: Raw data collected from the different sensors.
    - **Level 1**: Raw data stored to NetCDF file where attributes (such as sensors used, units, description of data, etc.) are added to the data. 
    Quality assurance is performed on the data and QA masks are generated.
    - **Level 2**: Processed data, this could include calculated parameters, transformed units, resampled or gridded data.

<font color='red'>*Add details about the data here.*</font>

<br />

### Folder `scripts`
Scripts used to read the data and perform further analysis.

<br />

### Folder `references`
Relevant scientific papers. 

<br />

### Folder `notes`
Anything to share (ideas, summaries, reports of analyses).

<br />

 ## Future updates

<font color='red'>*Add to-do list here*</font>

## Collaborators

- **Data collection**: Tomy Doda, Damien Bouffard, Bieito Fernández Castro,  Oscar Sepúlveda Steiner
- **Scripts and analysis**: Tomy Doda, Oscar Sepúlveda Steiner, Bieito Fernández Castro, Sebastiano Piccolroaz, Damien Bouffard

## Contact
- [Tomy Doda](mailto:tomy.doda@unil.ch)
- [Damien Bouffard](mailto:damien.bouffard@eawag.ch)


