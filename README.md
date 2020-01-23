<img src="https://github.com/geoportti/Logos/blob/master/geoportti_logo_300px.png">

# Spatiotemporal data fetching and processing scripts for 3-D lake modelling

This set of Python scripts automates the process of obtaining relevant geospatial
data for lake modelling purposes. This collection of scripts will
1) Fetch terrain, depth, weather, shoreline, river and hydrological data near 
the selected lake from openly available data sources.
2) Process terrain and depth data to a modelling grid.
3) Build weather and hydrological forcing input files.
4) Construct Fortran code for a basic lake model setup.

Some of the scripts can also be used independently for only accessing e.g.
weather or hydrological data.


## Dependencies and installation

### Puhti

To run the code on CSC's Puhti platform, install the code on location of your choosing. The
code has been tested with Python 3, but should also run under Python 2 with minor
changes.

The following extra environment might be needed:
```
module load geoconda
pip install xmljson --user
```

### Other platforms

The code has been tested on CSC's Puhti supercomputer, but should run on any
platform if required Python libraries and freely available datasets
are installed on the system.

At least the following extra Python libraries are needed:
```
dateutil  
geopandas  
numpy  
owslib  
pyproj  
pytz  
requests  
scipy  
shapely  
xmljson  
```


## Data sources

You will also need to make sure that some local datasets are installed.
On Puhti, most of the needed datasets are already installed for all users
under `/appl/data/geo`.


### Local datasets

1. A local copy of the map sheet information file, `utm25LR.shp`, is needed 
and can be downloaded from the
[NLS open data download service](https://www.maanmittauslaitos.fi/en/e-services/open-data-file-download-service).
This is also required on Puhti.

Following datasets are already installed on Puhti, but must be 
downloaded on other platforms:

2. The topographic database as map sheet zip files.
On Puhti, the database can be found under `/appl/data/geo/mml/maastotietokanta`.
On other platforms the user needs to download these from
[NLS open data download service](https://www.maanmittauslaitos.fi/en/e-services/open-data-file-download-service).

3. Lake shapefile (Ranta10 dataset) which on Puhti is under `/appl/data/geo/syke/ranta10jarvet`.
Can be downloaded from syke.fi --> Open information --> Spatial datasets --> Dataset packages --> Shoreline (lakes).

4. River network (Ranta10 dataset) which on Puhti is under `/appl/data/geo/syke/uomaverkosto`.
Can be downloaded from syke.fi --> Open information --> Spatial datasets --> Dataset packages --> Shoreline (river network).

### Remote datasets

The scripts fetch open data from following remote sources:

- TrafiCom's depth data for selected lakes from [Väylä's Web Feature Service](https://vayla.fi/web/en/open-data/download_services)
(WFS).

- Hydrological data from [Finnish Environment Institute's](https://www.syke.fi/en-US/Open_information/Open_web_services/Environmental_data_API)
(SYKE) OData service.

- Weather data from [Finnish Meteorological Institute's](https://en.ilmatieteenlaitos.fi/open-data)
(FMI) open WFS service.


## Configuration

The configuration file `config.py` is used when running the `gatherdata.py` script.
When running the code on Puhti, you only need to specify the path to the local 
map sheet file (`MAPSHEETFILE`). Other paths should contain working defaults.
In addition to paths and URIs, some other parameters can be fine-tuned. Further 
information can be found in the configuration file itself.

If you're running some of the scripts independently from `gatherdata.py`, you might
have to update some of the defaults inside those scripts.


## Scripts and files

Script                | Description
--------------------- | -----------
`build_model_code.py` | Auxiliary functions for building Fortran lake model code. Used by `gatherdata.py`.
`config.py`           | Configuration file for `gatherdata.py`.
`FMIwfs.py`           | Fetching and processing of weather data. Can be used independently.
`gatherdata.py`       | Main script for automatic lake model data gathering.
`ranta10.py`          | Functions for lake polygon processing. Used by `gatherdata.py`.
`SYKEdata.py`         | Fetching and processing of hydrological data. Can be used independently.
`VAYLAdata.py`        | Fetching and processing of depth data from Väylä. Can be used independently.

File                        | Description
--------------------------- | -----------
`fmi_stations_noAWOS.tsv`   | Tab-separated file with weather station metadata. Must be locally available. Used by `FMIwfs.py`. The file is generated from the web page at https://en.ilmatieteenlaitos.fi/observation-stations
`Usrdef_Model_template.f90` | COHERENS lake model Fortran code template. Used by `build_model_code.py`
`Usrdef_Surface_Data.f90`   | Example COHERENS lake model Fortran code for surface forcing.


## Lake model simulation code

Automatically generated [COHERENS](https://odnature.naturalsciences.be/coherens/)
lake model code works as a basis for a fully working 3-D lake model.
Full source code and manual for COHERENS can be downloaded from the project website.
These scripts only generate the `Usrdef_Model.f90` file which contains the basic
site specific model setup. *The actual lake model simulation environment
still needs to be set up by the user and is beyond the scope of this document*.

Inflow and outflow locations are automatically
detected and added to the model. However the vast majority of lakes do not have
all the necessary hydrological information available as measurements. Usually the
only available data is the outflow water level. Manual work is required to add
decent guesses of e.g. inflows.

Surface forcing (weather) is applied in the same way on every grid point.
This is a reasonable guess for most lakes. For the largest lakes and in some
use-cases a more detailed weather grid might be more desirable.


## Examples

It is assumed that the code and data are already installed and the necessary
environment loaded (see installation and configuration sections).

### Fetching and processing data for Lake Karhijärvi

1. Find out any coordinate within the lake in ETRS-TM35FIN coordinate system
from, for example, https://paikkatietoikkuna.fi/. In this case we use the
coordinates 6835478 N, 261628 E.

2. Edit `config.py` and set `POINTX = 261628.0` and `POINTY = 683478.0`.

3. Set `UTCTIME_MIN` and `UTCTIME_MAX` to the desired time interval. We use
`UTCTIME_MIN = "2016-07-01T00:00:00.000Z"` and 
`UTCTIME_MAX = "2016-08-01T00:00:00.000Z"` to fetch data from July, 2016.
The times must be specified in UTC.

4. Set `GRIDRESOLUTION` to the desired value in meters or use the default.
Using high resolutions makes the calculation much slower.

5. Feel free to play with the other settings too, but they are mostly for
fine-tuning the algorithms within the code.

6. Run the main script with `python gatherdata.py`

7. After the script has run, you should find the grid (`grid_coherens.dat`), 
weather data (`weather_coherens.dat`), hydrological data (`obdataXXX.dat`) 
and point cloud (`pointcloud.dat`) files saved on the current directory
as well as the `Usrdef_Model.f90` lake model setup file.

### Fetching weather data for Asikkala

1. Find out the desired coordinates in ETRS-TM35FIN.

2. Run the script `FMIwfs.py` with no parameters for help.

3. To get data nearest to the desired coordinates for the first two weeks
in July 2017 at 15 minute intervals:
`python FMIwfs.py 422025 6784263 2016-07-01T00:00:00+0000 2016-07-15T23:59:59+0000 15`

4. Data is fetched and possibly combined from multiple observation stations
and dumped to standard out (screen).

5. Note that the fetched data is always instantaneous (or nearly so) and is not averaged
to the fetched interval. It is possible to also fetch e.g. daily averages by slightly
modifying the code to ask for the correct values.


## Licensing

The data gathering and processing scripts are licensed with the MIT license.

The COHERENS Fortran code templates are licensed with the EUPL license.

Scripts written by Janne Ropponen, Finnish Environment Institute.


## Usage and Citing

When used, the following citing should be mentioned: "We made use of geospatial
data/instructions/computing resources provided by the Open Geospatial
Information Infrastructure for Research (oGIIR,
urn:nbn:fi:research-infras-2016072513) funded by the Academy of Finland."

