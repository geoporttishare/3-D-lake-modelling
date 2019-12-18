#!/usr/bin/python
#encoding: utf-8

from __future__ import print_function
from shapely.geometry import Point
import csv
import dateutil.parser
import numpy
import pytz
import requests
import shapely
import sys

##########################################################################
# Tools for downloading and processing basic hydrological data from SYKE's 
# Hydrology OData API for lake modelling purposes.
#
# https://www.syke.fi/fi-FI/content/37459/1#Hydrologiarajapinta
#
# oGIIR project
# Janne Ropponen/SYKE
# Last changed: 2019-12-17
##########################################################################

####################################################################
# getOData
####################################################################
def getOData(url, params=""):
   """Fetch data from an OData interface.

   Parameters:
   -----------
   url: string
        OData url
   params: string
        Query string

   Returns fetched data dictionary in json style.
   """
   r = requests.get(url, params=params)

   # Exit to system if fetching data fails in any way.
   if (r.status_code != requests.codes.ok):
      print("Cannot load data")
      print("URL: ", url)
      print("params: ", params)
      sys.exit(1)

   return r.json()

####################################################################
# getLocations
####################################################################
def getLocations(ODATA_URL_LOCATION, ETRSTM35FIN_XMIN, ETRSTM35FIN_XMAX, 
                                     ETRSTM35FIN_YMIN, ETRSTM35FIN_YMAX):
   """Fetch site location data from SYKE's OData API

   Parameters
   ----------
   ODATA_URL_LOCATION: string
        OData url for location
   ETRSTM35FIN_xxxx: float
        Bounding coordinates for locations

   Returns:
   List of location dictionaries that are within the query rectangle.
   E.g.
   [
    {
      "Paikka_Id": 1,
      "Suure_Id": 7,
      "Nimi": "J\u00e4nisjoki,Ruskeakoski ",
      "KoordErTmPohj":6926501,
      "KoordErTmIta":677699,
      ...
    }
   ]
   """
   # Build filter for locations inside the given coordinates
   filter  = "KoordErTmIta gt " + str(ETRSTM35FIN_XMIN) \
           + " and KoordErTmIta lt " + str(ETRSTM35FIN_XMAX) \
           + " and KoordErTmPohj gt " + str(ETRSTM35FIN_YMIN) \
           + " and KoordErTmPohj lt " + str(ETRSTM35FIN_YMAX) 
   payload = {'$filter' : filter}
   #
   data = getOData(ODATA_URL_LOCATION, payload)
   #
   locations = data['value']
   while 'odata.nextLink' in data.keys():
      data = getOData(data['odata.nextLink'])
      locations += data['value']
   return locations

####################################################################
# getvars
####################################################################
def getVars(ODATA_URL_VARS):
   """Fetch all offered variables from OData API url.

   Returns:
   List of variable metadata dictionaries. E.g.
   [
      {
         "Suure_Id": 1,
         "Suurekoodi": "W",
         "Nimi": "Vedenkorkeus",
         "Kuvaus": "Vedenkorkeus",
         "NimiEng": "Water level",
         "KuvausEng": "Water level",
         "Yksikko": "cm",
         "taulunnimi": "Vedenkorkeus",
         "DatataulunNimi": "Vedenkorkeus"
      }
   ]
   """
   data = getOData(ODATA_URL_VARS)
   vars = data['value']

   return vars

####################################################################
# generatelocationdatabase
####################################################################
def generateLocationDicts(locations, ODATA_URL_VARS, ODATA_ID_W, ODATA_ID_Q, 
                          searchpolygon):
   """Generates location dictionaries for water level (W) and
   discharge (Q). Dictionary key is location number (Paikka_Id) 
   and values are dictionaries with coordinates and location names.

   Parameters
   ----------
   locations: list of dictionaries
         Site locations and their metadata from OData API.
   ODATA_URL_VARS: 
         ODATA url for variable metadata.
   ODATA_ID_W: string
         OData water level name.
   ODATA_ID_Q: string
         OData discharge name.
   searchpolygon: polygon
         Area of interest.

   Returns W and Q location dictionaries.
   """
   varlist = getVars(ODATA_URL_VARS)
   locdb_W = {}
   locdb_Q = {}
   vardict = {}
   for var in varlist: # Build vartable where equivalent keys are the var id and it's code string
      vardict[var["Suure_Id"]] = var
      vardict[var["Suurekoodi"]] = var
   for loc in locations:
      d = {}
      d['x'] = loc['KoordErTmIta']
      d['y'] = loc['KoordErTmPohj']
      d['name'] = loc['Nimi']
      # Remove points too far outside lake polygon
      p = Point(d['x'],d['y'])
      if not p.intersects(searchpolygon):
         print("Skipped", d['name'], "(" + vardict[loc['Suure_Id']]['NimiEng'] + ") at coordinates", 
               d['x'], d['y'], "- outside search area.")
         continue
      # Distribute elevation and discharge data to separate location dicts
      if loc['Suure_Id'] == ODATA_ID_W: # Might want to change id lookup to var code lookup from vardict
         locdb_W[loc['Paikka_Id']] = d
         mtype = "elevation"
      elif loc['Suure_Id'] == ODATA_ID_Q:
         locdb_Q[loc['Paikka_Id']] = d
         mtype = "discharge"
      else:
         print("Skipped", d['name'], "(" + vardict[loc['Suure_Id']]['NimiEng'] + ") at coordinates",
               d['x'], d['y'])
         continue
      print("Found", mtype, "measurement site at", d['name'], "with coordinates", d['x'], d['y'])
   return locdb_W, locdb_Q

####################################################################
# processdata
####################################################################
def processData(dict, valuename, timename='Aika', multiplier=1.0):
   """Process "raw" OData dict and strip only the time and value.
   Also convert time to UTC and hydrodynamics model (COHERENS) format.

   Parameters
   ----------
   dict: dictionary
       Data dictionary as received from OData fetcher
   valuename: string
       Value nameto process
   timename: string
       Time field name
   multiplier: float
       Multiply value with this number. Useful in e.g. unit conversions.

   Returns dictionary with processed data.
   """
   # Gets valuename field from dict of sites along with timefield and multiplies values by multiplier
   # Returns dict of sites with list of values: time, coherenstime, value
   tz = pytz.timezone('Europe/Helsinki') # Default data timezone in case it doesn't exist
   if numpy.isnan(multiplier):
      print("Warning: multiplier ignored (NaN)")
      multiplier = 1.0
   newdict = {}
   for site in dict:
      newdata = []
      for meas in dict[site]:
         time = dateutil.parser.parse(meas[timename])
         # If timezone not present, assume local (Finland) timezone
         if time.tzinfo is None or time.tzinfo.utcoffset(time) is None:
            time = tz.localize(time)
         # If timezone is not UTC, convert time to UTC
         if time.tzname() != 'UTC':
            time = time.astimezone(pytz.utc)
         # Convert time from datetime object to COHERENS ASCII format
         coherenstime = time.strftime("%Y/%m/%d;%H:%M:%S,000")
         value = float(meas[valuename])*multiplier
         newdata.append([time, coherenstime, value])
      newdict[site] = newdata
   return newdict

####################################################################
# applyLevelCorr
####################################################################
def applyLevelCorr(Wdata, ODATA_URL, heightsystem):
   """Gets possible level correction data from OData API and
   applies it to each location Wdata.
   """
   for site in Wdata.keys():
      print("Applying", heightsystem, "level corrections to site id", site)
      filter  = "Paikka_Id eq " + str(site) \
              + " and TasoKoordinaatisto eq \'" + heightsystem + "\'" 
      payload = {'$filter' : filter}
      data = getOData(ODATA_URL, payload) # Example reply: {u'odata.metadata': u'http://rajapinnat.ymparisto.fi/api/Hydrologiarajapinta/1.1/odata/$metadata#VedenkTasoTieto', u'value': [{u'Tasokorjaus': 4408, u'TasoKoordinaatisto': u'N60', u'Korkeustaso_Id': 1, u'Paikka_Id': 2168}]}
      bias = float(data['value'][0]['Tasokorjaus']) 
      if bias != 0.0:
         for i, meas in enumerate(Wdata[site]):
            Wdata[site][i]['Arvo'] += bias
   return Wdata

####################################################################
# fetchOData
####################################################################
def fetchOData(locdb, ODATA_URL, UTCTIME_MIN, UTCTIME_MAX, text=""):
   """Fetches relevant data from SYKE's hydrological API

   Builds a proper request and fetches data from 
   ODATA_URL for sites in locdb and for time
   interval in UTCTIME_MIN--UTCTIME_MAX.
   Appends data from new links if the API says suggests
   that data continues.

   Parameters:
   locdb:       Location dictionary with keys as site ids.
   ODATA_URL:   OData request url (string)
   UTCTIME_MIN: Valid time string
   UTCTIME_MAX: Valid time string
   text:        Optional clarifying text string to display

   Returns a dictionary with site ids as keys and data as list.
   """
   datadict = {}
   trimmedlocdb = {} # Remove any sites that have no data
   for site in locdb.keys():
      print("Fetching "+text+"data for", locdb[site]['name'])
      filter  = "Paikka_Id eq " + str(site) \
              + " and Aika ge datetime\'" + UTCTIME_MIN + "\'" \
              + " and Aika le datetime\'" + UTCTIME_MAX + "\'"
      payload = {'$filter' : filter}
      # Get initial data based on "payload" request
      data = getOData(ODATA_URL, payload)
      values = data['value']
      # Get further data from possible new links
      while 'odata.nextLink' in data.keys():
         data = getOData(data['odata.nextLink'])
         values += data['value']
      if len(values)>0:
         datadict[site] = values
         print(" - found", len(values), "values for requested time interval.")
         trimmedlocdb[site] = locdb[site]
      else:
         print(" - no values for requested time interval.")
   return datadict, trimmedlocdb

####################################################################
# writeOData
####################################################################
def writeOData(locdict, datadict, prefix, unit):
   """Example code for writing fetched data to file in COHERENS
   hydrodynamic model compatible format.
   
   Parameters
   ----------
   locdict: dict
   datadict: dict
   prefix: string
       Filename prefix, e.g. "discharge"
   unit: string
       Data unit, e.g. "m3/s"
   """
   for site in datadict:
      filename = prefix + "_" + str(site) + ".dat"
      with open(filename,'w') as f:
         # Header (informative only)
         f.write( "#" 
                  + locdict[site]['name'] #.encode("utf-8") # encode needed if run in Python 2
                  + "; x=" 
                  + str(locdict[site]['x']) 
                  + ", y=" 
                  + str(locdict[site]['y']) 
                  + "\n" )
         f.write( "#timecode\t" + prefix + " [" + unit + "]\n")
         # Data
         for i in range(0,len(datadict[site])):
            f.write(str(datadict[site][i][1])+"\t"+str(datadict[site][i][2])+"\n" )
      print("Saved " + prefix + " data to "+filename)
   return

#############################################################################


####################################################################
# writeOBData
####################################################################
def writeOBData(locslist, W_data, Q_data, prefix="obdata", suffix=".dat"):
   """Write processed data to file in COHERENS hydrodynamic model 
   compatible format.
   
   Parameters
   ----------
   locslist: dict containing lists of u- and v-boundary data
   W_data, Q_data: dict
   prefix: string, optional
       Filename prefix, e.g. "obdata"
   suffix: string, optional
       Filename suffix, e.g. ".dat"
   """
   for i, site in enumerate(locslist['u']+locslist['v']):
      filename = prefix + str(i+1).zfill(3) + suffix
      with open(filename, 'w') as f:
         if site['id'] in W_data:
            type = 'W'
            data = W_data[site['id']]
         elif site['id'] in Q_data:
            type = 'Q'
            data = Q_data[site['id']]
         # Header (informative only, not needed by the model)
         f.write( "#" + site['name'] # .encode("utf-8") # encoded needed in Python 2
                  + "; x=" 
                  + str(site['x']) 
                  + ", y=" 
                  + str(site['y']) 
                  + "\n" 
                )
         f.write( "#timecode\t" + type + "_" + site['type'] + "\n" )
         # Data
         for j in range(0,len(data)):
            f.write(str(data[j][1])+"\t"+str(data[j][2])+"\n" )
      print("Saved " + type + "_" + site['type'] + " data to "+filename)
   return
# END writeOBData

#############################################################################

def getHydrologicalData(config, poly):
   """Wrapper for multiple functions to get and process OData discharge and elevation
   data within bounds.

   Parameters:
   config: dictionary
       Configuration object
   poly: shapely.Polygon
       Area where to search for data.
   """
   bounds = poly.bounds
   locations = getLocations(config.ODATA_URL_LOCATION, bounds[0], bounds[2], bounds[1], bounds[3])

   W_locs, Q_locs = generateLocationDicts(locations, 
                                          config.ODATA_URL_VARS, 
                                          config.ODATA_ID_W, 
                                          config.ODATA_ID_Q, 
                                          poly)

   # Get and process water stage data
   W_rawdata, W_locs = fetchOData(W_locs, config.ODATA_URL_W, config.UTCTIME_MIN, config.UTCTIME_MAX, "water stage ")
   W_rawdata = applyLevelCorr(W_rawdata, config.ODATA_URL_W_LEVELCORR, config.HEIGHTSYSTEM)
   W_data = processData(W_rawdata, 'Arvo', multiplier = 0.01) # Change units from cm to m

   # Get and process discharge data
   Q_rawdata, Q_locs = fetchOData(Q_locs, config.ODATA_URL_Q, config.UTCTIME_MIN, config.UTCTIME_MAX, "discharge ")
   Q_data = processData(Q_rawdata, 'Arvo') # Process raw data to usable format

   return W_locs, Q_locs, W_data, Q_data
# END getHydrologicalData
   
#############################################################################

########
# MAIN #
########
def main():
   """Demonstration code to get waterlevel and discharge data for
   requested area and time interval. Saves data in COHERENS model
   format. Configuration and request parameters below.
   """
   if len(sys.argv)!=7:
      print("Fetch hydrological data SYKE OData service.")
      print("Usage: python "+sys.argv[0]+" <xmin> <ymin> <xmax> <ymax> <timeUTCstart> <timeUTCend>")
      print("Example (Oulujarvi region): python "+sys.argv[0]+" 482900 7107700 554000 7160000 2016-07-09T00:00:00Z 2016-07-15T23:59:59Z")
      print("xmin etc are bounding box coordinates in ETRS-TM35FIN.")
      sys.exit(1)

   # CONFIGURATION
   ODATA_BASE_URL = "http://rajapinnat.ymparisto.fi/api/Hydrologiarajapinta/1.1/odata/"
   ODATA_URL_LOCATION = ODATA_BASE_URL + "Paikka" 
   ODATA_URL_W = ODATA_BASE_URL + "Vedenkorkeus"
   ODATA_URL_W_LEVELCORR = ODATA_BASE_URL + "VedenkTasoTieto"
   ODATA_URL_Q = ODATA_BASE_URL + "Virtaama"
   ODATA_URL_VARS = ODATA_BASE_URL + "Suure"
   ODATA_ID_W = 1 # Suure_id in OData for elevation
   ODATA_ID_Q = 2 # Suure_id in OData for discharge
   #HEIGHTSYSTEM = "N60" # Get data in this elevation system

   # Coordinate rectangle for Pyhajarvi (Sakyla), ETRS-TM35FIN
   ETRSTM35FIN_XMIN = float(sys.argv[1])
   ETRSTM35FIN_YMIN = float(sys.argv[2])
   ETRSTM35FIN_XMAX = float(sys.argv[3])
   ETRSTM35FIN_YMAX = float(sys.argv[4])

   # Request time interval
   UTCTIME_MIN = sys.argv[5]
   UTCTIME_MAX = sys.argv[6]

   # Get all measurement locations that match the coordinates and time
   locations = getLocations(ODATA_URL_LOCATION, ETRSTM35FIN_XMIN, ETRSTM35FIN_XMAX, 
                            ETRSTM35FIN_YMIN, ETRSTM35FIN_YMAX)

   # Get data locations for water stage W and discharge Q
   poly = shapely.geometry.box(ETRSTM35FIN_XMIN, ETRSTM35FIN_YMIN, ETRSTM35FIN_XMAX, ETRSTM35FIN_YMAX)
   locdb_W, locdb_Q = generateLocationDicts(locations, ODATA_URL_VARS, ODATA_ID_W, ODATA_ID_Q, poly)

   # Get W, Q measurement data from ODATA API and store it in a dictionary (key = site)
   W, locdb_W = fetchOData(locdb_W, ODATA_URL_W, UTCTIME_MIN, UTCTIME_MAX, "water stage ")
   Q, locdb_Q = fetchOData(locdb_Q, ODATA_URL_Q, UTCTIME_MIN, UTCTIME_MAX, "discharge ")

   # Clean fetched data and process to usable format
   newQ = processData(Q, 'Arvo', multiplier = 1.0)
   newW = processData(W, 'Arvo', multiplier = 0.01)

   # Save data to file
   writeOData(locdb_Q, newQ, "discharge", "m3/s")
   writeOData(locdb_W, newW, "waterstage", "m")


# Run from shell
if __name__ == "__main__":
   main()

