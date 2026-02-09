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

### 1. Matlab installation

Matlab is required to run some of the scripts.

### 2. Python installation

Python 3 is required to run some of the scripts. Three types of installation are possible:
- Recommended option: download [Miniforge](https://github.com/conda-forge/miniforge). 
- User-friendly option: download the [Anaconda distribution](https://www.anaconda.com/products/individual).
- Classic option: download Python from the [official website](https://www.python.org/downloads/).

### 3. Repository installation

- If using GIT, clone the repository to your local machine using the command in Git Bash: 

    ``` 
    git clone https://github.com/tdoda/turbulence_methods.git
    ```
 
    Note that the repository will be copied to your current working directory.
- Without GIT, just download the entire ZIP folder from https://github.com/tdoda/turbulence_methods.git ("Code" > "Download ZIP") and extract it.

### 4. Python packages installation

<font color='red'>*Need to update environment and requirements files.*
</font>

1. Open the terminal (e.g., Anaconda Prompt), and move to the `turbulence_methods` repository.
2. Create a new environment *turbulence-methods* and install the packages as follows:
    - If using conda (Anaconda or Miniforge installation):
        ```
        conda env create -f environment.yml
        conda activate turbulence-methods 
        ```
        It is also possible to install the packages from `requirements.txt` with pip instead:
        ```
        conda create -n turbulence-methods python=3.11
        conda activate turbulence-methods
        pip install -r requirements.txt
        ```
    - If using mamba (Anaconda or Miniforge installation):
        ```
        mamba env create -f environment.yml
        mamba activate turbulence-methods 
        ```
    - If using pip (classic Python installation):
        ```
        python -m venv turbulence-methods       
        source turbulence-methods /bin/activate  # For Linux/macOS
        turbulence-methods\Scripts\activate     # For Windows
        pip install -r requirements.txt
        ```

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

## Example data

<font color='red'>*Add information about the example datasets (geospatial location, temporal coverage).*</font>

## Usage
### Processing data from microstructure profilers
#### Pre-deployment setup
1. Before any deployment, create a `Level0` folder on the field laptop, where the data of the field campaign will be saved.
2. Save the configuration file (`*.cfg`) used by the profiler in the `Level0` folder:
    - VMP: copy the configuration file from a previous field campaign.
    - microCTD: extract the configuration file from the profiler using RSILink (see microCTD SOP).

    It is possible to modify the configuration file before deployment by specifying the correct probes serial numbers and shear sentivities, but it is not mandatory as everything can be done post-deployment. If the configuration file is modified before deployment, **you need to load it to the profiler after modification**. This is anyway done via ODAS5RT for the VMP, but it has to be done with RSILink for the microCTD (see SOP). The configuration file stored in the `Level0` folder **must always** correspond to the configuration loaded to the profiler during the campaign. Note that such configuration file can anyway be retrieved later from the `*.P` files, but it is a good practice to save it separately.
3. Prepare the file `scripts/microstructure/load_parameters_*lake_name*.m` (replacing *lake_name* by the name of the lake) by creating a new section for the field campaign. Parameters that can already be written are `param.folder`, `param.instrument`, `param.prof_dir` and `param.config` (see `notes/Turbulence_methods_manual.pdf` for more details on the parameters).

    For example:

    ```
    elseif strcmp(date,'20251127') 
        param.folder = [general_data_folder,'20251127\Level0\'];
        param.filename_list={}; % Several files can be listed here

        % Set P offset and sh probe sensitivity
        % param.offset_P=;
        % Use shear sensitivities specified in config file: 
        param.cfgfile = '';

        param.instrument='VMP';
        param.info.prof_dir = 'down';
        param.config.T1=true;
        param.config.T2=true;
        param.config.S1=true;
        param.config.S2=true;
        param.config.uC1=false;
        param.config.uC2=false;
    ```

4. Prepare the file `scripts/analyze_VMP_microCTD/analyze_profiles.m` by specifying the lake name, the campaign name matching the new section in `scripts/microstructure/load_parameters_*lake_name*.m` and the paths to the different folders:
    ```
    lakename='Zug'; 
    
    general_data_folder='..\..\data\microCTD\'; % Where fieldwork data is stored
    
    odas_folder='..\..\functions\odas_v4.4\'; % Where ODAS functions are stored
    
    functions_folder="..\..\functions\microstructure\"; % Where microstructure functions are stored

    date_campaign="20251127"; 
    ```
Make also sure that the correct `load_parameters_*lake_name*.m` function is called:
```
param=load_parameters_Zug(lakename,date_campaign,general_data_folder);
```


#### Deployment: collecting and saving data

Make sure you can have access to Matlab in the field (i.e., VPN connection).

Save the data in the `Level0` folder: 
- VMP: when starting the instrument, select the `*.cfg` file stored in the `Level0` folder via ODAS5RT. Data will be automatically saved in this folder.  
- microCTD: after deployment, download the data and save it to the `Level0` folder (see microCTD SOP).

#### Post-deployment parameters

In `scripts/microstructure/load_parameters_*lake_name*.m`, update the following parameters required to import and process the data (following two sections):
- *filename_list* (mandatory): name of Level0 files to process without the file extension. For example:
    ```
    param.filename_list={'DAT_053','DAT_055'};
    ```
- *cfgfile* (optional): specify the configuration file used for data processing (without the file extension), containing the correct shear sensitivities and pressure offset (if the "offset" method is used for pressure correction). The easiest option is to copy the original configuration file in the `Level0` folder, modify the required parameters and save it as an updated configuration file. For example:
    ```
     param.cfgfile = 'SETUP_updated';
    ```
    If *cfgfile* is left empty, the shear sensitivities **must** be specified as parameters:
    ```
    # Example:
    param.S_sh1=0.0603 # [V.s^2.m^(-2)]
    param.S_sh2=0.0824 # [V.s^2.m^(-2)]
    ```
    
    except if the configuration of the deployment is used for data processing by specifying
    ```
    modif_cfg=false
    ```
    in `scripts/analyze_VMP_microCTD/analyze_profiles.m`. 
    
    If the "offset" method is used for pressure correction, the pressure offset must either be included in the coef0 of the pressure channel in the configuration file (if *cfgfile* is provided), or specified as a parameter (if *cfgfile* is empty):
    ```
    param.offset_P=0.1; # [dbar]
    ```
    or entered manually while running `scripts/analyze_VMP_microCTD/analyze_profiles.m` (if *cfgfile* and *P_offset* are empty).

- *atm_press_method* (optional): pressure correction method, with one of the following values:
    - *offset*: apply a pressure offset from the configuration file or with the parameter *offset_P*, typical method for VMP
    - *min*: minimum air pressure in the data file is considered as the atmospheric pressure, works only if each data files starts in the air (typical for microCTD)
    - *cond*: use microconductivity to detect the surface, works only for upward profiles
    - *FP07*: use FP07 to detect the surface, works only for upward profiles.
    
    If *atm_press_method* is empty or inexistant, the "offset" method will be used by default.

#### Checking the profiles 

To have a quick look at the collected profiles **without** performing the turbulence analysis, run `scripts/analyze_VMP_microCTD/analyze_profiles.m` with
```
turbulence_analysis=false;
```
and specify whether the data and figures should be saved:
```
save_checkdata=true; % If =true, save the "checked" data
save_checkfig=true; % If =true, save the "checked" figures
```


Other options to specify are:
- *modify_cfg*: whether the configuration file should be modified with the correct shear sensitivities → required if the shear sentivities were not included in the pre-deployment configuration file:
    ```
    modify_cfg=true;
    ```
- *calibrate_FP07*: whether the FP07 data should be calibrated based on the CTD data → required :
    ```
    calibrate_FP07=true;
    ```

#### Running the turbulence analysis

To perform the full turbulence analysis, run `scripts/analyze_VMP_microCTD/analyze_profiles.m` with
```
turbulence_analysis=true;
```

Data (`*.netCDF` and `*.mat` files) and figures will be exported to Level1 and Level2 folders.

Other options to specify are:
- *modify_cfg* (usually *true*): whether the configuration file should be modified with the correct shear sensitivities → required if the shear sentivities were not included in the pre-deployment configuration file.
    ```
    modify_cfg=true;
    ```
- *calibrate_FP07* (usually *true*): whether the FP07 data should be calibrated based on the CTD data → required.
    ```
    calibrate_FP07=true;
    ```
- *run_quick_look* (usually *false*): whether the `quick_look()` ODAS function should be run to perform shear turbulence analysis from Rockland. This could be a way of comparing the shear-based turbulence estimates.
    ```
    run_quick_look=false;
    ```
- *run_dissip* (usually *true*): whether the shear and temperature-based turbulence analysis from Piccolroaz et al. (2021) should be performed.
    ```
    run_dissip=true;
    ```
- *make_plot_prof * (usually *true*): whether profile plots should be made.
    ```
    make_plot_prof=true;
    ```
- *ind_plot_spectra* (usually empty): indices of bins where spectrum should be plotted.
    ```
    # Example:
    ind_plot_spectra = [5,10];
    ```
- *show_progress* (usually *true*): whether the progress should be displayed in the terminal.
    ```
    show_progress=true;
    ```

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


