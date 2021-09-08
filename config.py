################################################################################################
# CONFIGURATION
################################################################################################

# Paths set up for running in Puhti.csc.fi

######################
# LOCAL data sources
######################

# Lake shapes (Ranta10 Järvet). Can be downloaded from syke.fi --> Avoin tieto.
SYKERANTA10SHP = "/appl/data/geo/syke/ranta10jarvet/jarvi10.shp"

# Map sheet divisions. Can be downloaded from MML open data file service.
MAPSHEETFILE = "/users/jroppone/ogiir/karttaruudut/utm25LR.shp"

# Terrain database. Can be downloaded from MML open data file service. 
# Zipped shapefiles are saved in a directory structure based on the map sheet names.
# E.g. TERRAINDBDIR/N4/N42/N4211L.shp.zip
#      TERRAINDBDIR/N4/N42/N4211R.shp.zip 
#      etc.
TERRAINDBDIR = "/appl/data/geo/mml/maastotietokanta/2018"

# Waterway/channel shapes (Ranta10 Uomaverkosto). Can be downloaded from syke.fi --> Avoin tieto.
CHANNELDB = "/appl/data/geo/syke/uomaverkosto/Uoma10.shp"

####################################
# REMOTE data sources and settings
####################################

# TrafiCom depth data WFS URL
TRAFICOMURL = "https://julkinen.traficom.fi/inspirepalvelu/rajoitettu/wfs?request=getcapabilities"

# SYKE hydrological database OData URL and settings
ODATA_BASE_URL = "http://rajapinnat.ymparisto.fi/api/Hydrologiarajapinta/1.1/odata/"
ODATA_URL_LOCATION = ODATA_BASE_URL + "Paikka" 
ODATA_URL_W = ODATA_BASE_URL + "Vedenkorkeus"
ODATA_URL_W_LEVELCORR = ODATA_BASE_URL + "VedenkTasoTieto"
ODATA_URL_Q = ODATA_BASE_URL + "Virtaama"
ODATA_URL_VARS = ODATA_BASE_URL + "Suure"
ODATA_ID_W = 1 # Suure_id in OData for elevation
ODATA_ID_Q = 2 # Suure_id in OData for discharge

# FMI weather data WFS URL and settings
FMIWFSURL = 'https://opendata.fmi.fi/wfs'
FMIPARAMS = ['temperature', 'windspeedms', 'winddirection',
             'humidity', 'pressure', 'ri_10min', 'totalcloudcover']
FMISTATIONSFILE = 'fmi_stations_noAWOS.tsv'
FMIMAXREQLENGTH = 7 # Request max 7 days of data at a time
FMITIMESTEP = 10 # (minutes) Min & Max time must be consistent with timestep, i.e. divisible

######################
# GENERAL SETTINGS
######################

UTCTIME_MIN = "2021-07-01T00:00:00.000Z" # Temporal data fetch start time
UTCTIME_MAX = "2021-08-01T00:00:00.000Z" # Temporal data fetch end time

POINTCLOUDFILE = "pointcloud.dat" # Output file for point cloud

HEIGHTSYSTEM = "N60" # Get data in this elevation system (terrain database uses N60)

AREATOLERANCE = 1.0 # [m2] Tolerance in area calculations
GRIDRESOLUTION = 50.0 # [m] Depth grid resolution
BUFFERZONE = 100.0 + 3.0*GRIDRESOLUTION # Bufferzone outside the lake within which height data is collected. Also acts as a guarantee that we get correct MWL and lake name points within the loaded map sheets.
SEARCHRADIUS = 500.0 # [m] Radius when searching for MWL and name points for lake polygon
CHANNELSEARCHRADIUS = 10.0 # [m] Radius from lake polygon border when searching for incoming/outgoing rivers
INTERPOLATIONMETHOD = 'linear' # Grid generation interpolation method. Other options: cubic, nearest

MAXCOMPLEXITY = 1000 # Maximum lake complexity before simplifying geometry (for performance purposes)
SIMPLIFYTOLERANCE = 100.0 # [m] Simplifying tolerance.

#####################################################
# COHERENS lake simulation code generation settings
#####################################################
COHERENS_MAXOUTTIME = 600 # [s] Make sure output is possible at least this often.
COHERENS_IC3D = 5 # 3-D timestep
COHERENS_MAXDELT2D = COHERENS_MAXOUTTIME/COHERENS_IC3D # MAXDELT2D must be divisible by 2d-timestep
COHERENS_EPSILON = 10E-6
COHERENS_NZ = 10 # Sigma layers
COHERENS_WEATHEROUTFILE = "weather_coherens.dat"
COHERENS_GRIDFILE = "grid_coherens.dat"
COHERENS_USRDEFMODEL_TEMPLATE = 'Usrdef_Model_template.f90'
COHERENS_USRDEFMODEL = 'Usrdef_Model.f90'
   

######################
# Default coordinates
######################

#POINTX = 334068.0 # Lappajarvi
#POINTY = 7008459.0

# Other example lake coordinates

#POINTX = 500292.0 # Iisvesi x
#POINTY = 6948679.0 # Iisvesi y

#POINTX = 482668.0 # Puula
#POINTY = 6847723.0 # Puula

#POINTX = 261559.0 # Karhijarvi
#POINTY = 6835970.0 # Karhijarvi

#POINTX = 543122.0 # Syvari (only LiVi depths)
#POINTY = 7032548.0

#POINTX = 509056.0 # Oulujarvi
#POINTY = 7129280.0

#POINTX = 552178.0 # Nuasjarvi
#POINTY = 7116953.0

#POINTX = 552109.0 # Unnukka
#POINTY = 6915599.0

#POINTX = 520510.0 # Kemijarvi
#POINTY = 7400988.0 #

#POINTX = 511853.0 # Heinanen
#POINTY = 6950931.0 # Heinanen

POINTX = 246111.0 # Sakylan Pyhajarvi
POINTY = 6772561.0

#POINTX = 436446.0 # Jyvasjarvi
#POINTY = 6901246.0

#POINTX = 284691.0 # Rautavesi
#POINTY = 6808657.0

#POINTX = 302396.0 # Kulovesi
#POINTY = 6819888.0

# POINTX = 377144.0 # Pyhajarvi (Hauho)
# POINTY = 6783843.0

