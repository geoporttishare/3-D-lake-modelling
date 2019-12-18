#!/usr/bin/python
from __future__ import print_function
from dateutil import parser as dateparser
from shapely.geometry import Point
from owslib.wfs import WebFeatureService
import copy
import csv
import dateutil
import datetime
import sys
import json
import math
import numpy
import pyproj
import re
import xmljson
import xml.etree.ElementTree

##########################################################################
# Downloads FMI weather data from their open Web Feature Service (WFS)
#
# Designed for getting data for lake models.
#
# Some features:
#  * Finds closest measurement station to given coordinates
#  * If data is not available for wanted time period on given station, 
#    gets data from next closest station
#  * Combines data from different stations
#  * Basic data interpolation for NaN values (not perfect)
#  * Writes data in COHERENS model format
#
# oGIIR project
# Janne Ropponen/SYKE
# Last changed: 2019-12-12
##########################################################################

# Improvements TO DO:
# * Fill (interpolate) gaps
# * Postprocess msl atmpres to local atmpres, wind to components, rain to mm/s etc.
# * Rain NaN = 0.0

def dttostr(dt):
   """Converts datetime to string
   Input: datetime object
   Output: date string
   """
   return dt.strftime("%Y-%m-%dT%H:%M:%S%z")

def strtodt(string):
   """Converts string to datatime object
   Input:  date string 
   Output: datetime object
   """
   return dateparser.parse(string)

def readWeatherStationsFile(stationsfile, x, y):
   """Reads a local copy of FMI station list and retains 
   only weather and precipitation stations.

   Returns a stations dictionary with FMISIDs as keys.
   """
   stations = {}
   proj_etrstm35fin = pyproj.Proj(init='epsg:3067')
   proj_wgs84 = pyproj.Proj(init='epsg:4326')
   # Need to project wgs84 coordinates to etrs-tm35fin in order to reliably compare distances
   # because shapely distance expresses the Euclidean distance between points.
   p = Point(x,y) # Lake location in etrs-tm35fin
   # TODO: Get x,y from the middle of the lake/shape, not the place user specified.
   with open(stationsfile) as stfile:
      contents = csv.reader(stfile, delimiter='\t')
      next(contents) # skip header
      for row in contents:
         types = str.lower(row[7]).split(';')
         if 'weather' in types or 'precipitation' in types:
            stations[row[1]] = {}
            stations[row[1]]['name'] = row[0]
            stations[row[1]]['latitude'] = row[4]
            stations[row[1]]['longitude'] = row[5]
            stations[row[1]]['startyear'] = row[8]
            stations[row[1]]['endyear'] = row[9]
            stations[row[1]]['types'] = types
            coords = pyproj.transform(proj_wgs84,proj_etrstm35fin,row[5],row[4])
            distance = p.distance(Point(coords))
            stations[row[1]]['xloc'] = coords[0]
            stations[row[1]]['yloc'] = coords[1]
            stations[row[1]]['distance'] = distance # Distance to Point p in meters
   print("Found", len(stations), "weather/precipitation stations in", stationsfile)
   return stations

def filterStations(stations, timeStart, timeEnd):
   """Returns stations filtered by operation dates and sorted by distance.
   """
   filteredstations = {}
   for st in stations:
      year0 = int(stations[st]['startyear'])
      yearnstr = stations[st]['endyear']
      if yearnstr == '':
         yearn = timeEnd.year
      else:
         yearn = int(yearnstr)
      # Return only stations that might have data on the time interval requested
      if year0<=timeEnd.year and yearn>=timeStart.year:
         filteredstations[st] = stations[st]
   # Make sorted dictionary based on distance
   distlist = []
   for fmisid in filteredstations:
      distlist.append( [fmisid, filteredstations[fmisid]['distance']] )
   distlist.sort(key=lambda x: x[1]) # sort by distance
   # Add distance number (closest to farthest) to stations dictionary
   sortedstations = {}
   for i, st in enumerate(distlist,1):
      sortedstations[i] = {}
      sortedstations[i]['fmisid'] = st[0]
      for key in filteredstations[st[0]]:
         sortedstations[i][key] = filteredstations[st[0]][key]
      filteredstations[st[0]]['closest'] = i # Closest is number 1, farthest is number n
   print("Found", len(filteredstations), 
         "stations with data on requested interval. Closest is",
         sortedstations[1]['distance']/1000.0, "km from target area.")
   return filteredstations, sortedstations

def fmiwfsrequest(fmiwfs, params, timestep, fmisid, timeStart, timeEnd, maxreqlength):
   """
   Requests data from FMI WFS service in manageable chunks and processes
   the fetched data to a dictionary.
   """
   times = []
   # Split requests to maxreqlength intervals if needed.
   if (timeEnd-timeStart)>datetime.timedelta(days=maxreqlength):
      totsecs = (timeEnd-timeStart).total_seconds()
      maxreqsecs = datetime.timedelta(days=maxreqlength).total_seconds()
      nbintervals = int(math.ceil(totsecs/maxreqsecs))
      for n in range(0,nbintervals):
         time0 = timeStart + datetime.timedelta(days=n*maxreqlength)
         time1 = timeStart + datetime.timedelta(days=(n+1)*maxreqlength)
         times.append( (time0, min(time1,timeEnd)) )
   else:
      times.append((timeStart,timeEnd))
   # Get and process data from FMI WFS service
   alltimes = numpy.empty((0,0))
   alldata = numpy.empty((0,len(params))) # second dimension size == number of parameters
   paramstr = ','.join(params)
   print("Accessing FMI WFS for data from", times[0][0], "to", times[-1][1])
   for t, interval in enumerate(times,1):
      response = fmiwfs.getfeature( 
                    storedQueryID = 'fmi::observations::weather::multipointcoverage', 
                    storedQueryParams = {
                       'parameters' : paramstr,
                       'fmisid'     : fmisid,
                       'starttime'  : dttostr(interval[0]), 
                       'endtime'    : dttostr(interval[1]),
                       'timestep'   : timestep 
                    }
                 )
      print(round(float(t)/len(times)*100),"% of data fetched.")
      #
      nbvars, varnames, fmisid_wfs, station, latlon, newtimes, newdata = processWFSResponse(response) # Get actual data from xml response
      if nbvars is not None:
         # Add this request's data to full array
         alltimes = numpy.append(alltimes,newtimes)
         alldata = numpy.append(alldata,newdata,axis=0)
   # Reorganize data to a dictionary
   datadict = {}
   datadict['fmisid'] = fmisid
   datadict['station'] = station
   datadict['nbvars'] = nbvars
   datadict['varnames'] = varnames
   datadict['latlon'] = latlon
   datadict['length'] = len(alltimes)
   datadict['times'] = alltimes
   datadict['values'] = alldata
   #
   return datadict

def processWFSResponse(response):
   """Used to convert XML response to relevant lists of data.
   Expects XML response in MultiPointCoverage format.
   Yes, it's ugly.
   """
   # For reasons mysterious, sometimes the response is not 
   # a cStringIO object but a ResponseWrapper object.
   try:
      xmlstring = response.getvalue()
   except AttributeError: 
      xmlstring = response.read()
   # A hack to remove namespaces and gml attributes from XML response
   # to make further processing easier. Yeah, not proud about this. But it works.
   xmlstring = re.sub('<[a-zA-Z0-9]*?:', '<',  xmlstring) # Strip namespace:xxx attributes from all tag beginnings
   xmlstring = re.sub('</.*?:',          '</', xmlstring) # Strip namespace:xxx attributes from all tag ends
   xmlstring = re.sub('gml:',            '',   xmlstring) # Strip gml: attributes from all tag ends
   #
   # Convert data from xml string to json to dictionary using the BadgerFish notation
   bf = xmljson.BadgerFish()
   data = bf.data(xml.etree.ElementTree.fromstring(xmlstring))
   #
   # Extract relevant data from data dictionary
   if int(data['FeatureCollection']['@numberReturned']) > 0:
      fmisid = data['FeatureCollection']['member']['GridSeriesObservation'] \
                   ['featureOfInterest']['SF_SpatialSamplingFeature'] \
                   ['sampledFeature']['LocationCollection']['member'] \
                   ['Location']['identifier']['$']
      #
      station = data['FeatureCollection']['member']['GridSeriesObservation'] \
                    ['featureOfInterest']['SF_SpatialSamplingFeature'] \
                    ['sampledFeature']['LocationCollection']['member'] \
                    ['Location']['name'][0]['$']
      #
      latlontime = data['FeatureCollection']['member']['GridSeriesObservation'] \
                       ['result']['MultiPointCoverage']['domainSet'] \
                       ['SimpleMultiPoint']['positions']['$']
      #
      values = data['FeatureCollection']['member']['GridSeriesObservation'] \
                   ['result']['MultiPointCoverage']['rangeSet']['DataBlock'] \
                   ['doubleOrNilReasonTupleList']['$']
      #
      # Get names of variables
      varlist = data['FeatureCollection']['member']['GridSeriesObservation'] \
                    ['result']['MultiPointCoverage']['rangeType']['DataRecord'] \
                    ['field']
      #
      if len(varlist[0])==1: # Try to recognise if returned data is STRING instead of a LIST
         varlist = [varlist] # Convert string to list with only element being the string
      nbvars = len(varlist)
      varnames = []
      for var in varlist:
         varnames.append(var['@name'])
      #
      # Convert data string to numpy array and reshape 
      valuearray = numpy.fromstring(values,sep=' ') # parse data to 1-D array
      valuearray = valuearray.reshape((len(valuearray)//nbvars,nbvars)) # Reshape 1-D array to 2-D array
      coordarray = numpy.fromstring(latlontime,sep=' ') # parse data to 1-D array
      coordarray = coordarray.reshape((len(coordarray)//3,3)) # Reshape 1-D array to 2-D array
      times = coordarray[:,2]
   else:
      print("Zero matches returned")
      return None, None, None, None, None, None, None
   #
   return nbvars, varnames, fmisid, station, coordarray[0,0:2], times, valuearray

def mergeFMIdata(basedict, newdict):
   """Merge data from two FMI data dictionaries. Assumes basedict 
   has all the possible vars and newdict replaces some of them.
   Also assumes that 'times' are exactly equal in both dictionaries.

   TO DO: save which varnames belong to which station.
   """
   # Append new station data to basedict data
   basedict['latlon'].append(newdict['latlon'])
   basedict['fmisid'].append(newdict['fmisid'])
   basedict['station'].append(newdict['station'])
   # Copy columns from newdict to replace the same columns in basedict
   for n, newvar in enumerate(newdict['varnames']):
      for b, basevar in enumerate(basedict['varnames']):
         if basevar==newvar:
            basedict['values'][:,b] = newdict['values'][:,n]
   return basedict

def cleanFMIdata(fmidata):
   """Clean up NaNs from data and interpolate/extrapolate when necessary.
   """
   print("Interpolating weather data.")
   # 1. Precipitation: Replace NaNs with zeroes in-place
   try:
      rainind = fmidata['varnames'].index('ri_10min') # unit = mm/h
      tmp = numpy.nan_to_num(fmidata['values'][:,rainind],copy=False)
   except ValueError:
      print("Warning: Rain column (ri_10min) not in array.")
      pass
   # 2. Interpolate NaNs out of the array
   for v in range(0,len(fmidata['varnames'])):
      notnans = ~numpy.isnan(fmidata['values'][:,v])
      x = fmidata['times']
      xp = fmidata['times'][notnans]
      fp = fmidata['values'][:,v][notnans]
      # Prepare column for circular interpolation in case of wind direction
      # https://stackoverflow.com/questions/27295494/bounded-circular-interpolation-in-python
      # deg2rad + unwrap + rad2deg + %360
      if fmidata['varnames'][v] == 'winddirection':
         fp = numpy.rad2deg(
               numpy.unwrap(
                numpy.deg2rad(fp)
               )
              )
         fmidata['values'][:,v] = numpy.interp(x,xp,fp)%360
      else:
         fmidata['values'][:,v] = numpy.interp(x,xp,fp)
   #
   return fmidata
   # 3. To do:
   #  - cloudiness to fraction, relative humidity to fraction
   
def getFMIdata(FMIWFSURL, FMIPARAMS, FMISTATIONSFILE, MAXREQLENGTH,
               UTCTIME_MIN, UTCTIME_MAX, TIMESTEP,
               POINTX, POINTY):
   """Returns weather data nearest to user supplied coordinates 
   for user supplied time period.
   """
   # Operation:
   # 1) Get data from closest station that is in operation during the requested time interval
   # 2) If some column of data is missing (nan), try to get same data from the next nearest 
   # station in operation during the requested time interval
   # 3) Continue fetching missing data from stations farther away until we get it all, no more 
   # data is available, or we hit some predefined distance limit (e.g. 100 km)
   # 4) Save metadata about fetched data locations
   # TO DO: Make a script to update local fmi_stations.tsv automatically
   #
   timeStart = strtodt(UTCTIME_MIN)
   timeEnd = strtodt(UTCTIME_MAX)

   # Build a weather stations dictionary from file
   stations = readWeatherStationsFile(FMISTATIONSFILE,POINTX,POINTY)

   # Get only relevant stations and sort them by distance
   filteredstations, sortedstations = filterStations(stations, timeStart, timeEnd)
   fmiwfs = WebFeatureService(url=FMIWFSURL, version='2.0.0')

   # Fetch data from FMI WFS and save the results to a stationwise list of data dictionaries
   fulldata = [] # A list of dictionaries containing fetched measurement data
   fmiparams = copy.deepcopy(FMIPARAMS) # Make a copy of parameter list  because the list will be modified
   for i in range(1,len(sortedstations)+1):
      if len(fmiparams)==0: # Data for all parameters acquired, stop going through stations
         break
      print("Requesting data from station",sortedstations[i]['name'],"fmisid =",sortedstations[i]['fmisid'])
      datadict = fmiwfsrequest(fmiwfs, fmiparams, TIMESTEP, sortedstations[i]['fmisid'], timeStart, timeEnd, MAXREQLENGTH)
      if len(datadict['values'])==0: # No data returned at all
         print(" - no data returned for requested time interval, retrying from next closest station.")
         continue
      else:
         # At least some data was returned. Now check that all variables have at least some non-nan data
         # and if not, try to get the missing variable data from the next closest fmi station
         newparams = []
         for j in range(0,len(fmiparams)):
            if numpy.isnan(datadict['values'][:,j]).all():
               newparams.append(fmiparams[j])
               print(" - no data for parameter",datadict['varnames'][j],"- retrying from next closest location.")
         fmiparams = newparams
      fulldata.append(datadict)

   # Information
   print("Data fetched from",len(fulldata),"different stations.")
   for i in range(0,len(fulldata)):
      print("Total of",len(fulldata[i]['values']),"timesteps fetched from station",fulldata[i]['station'])
   fulldata[0]['latlon'] = [fulldata[0]['latlon']]
   fulldata[0]['fmisid'] = [fulldata[0]['fmisid']]
   fulldata[0]['station'] = [fulldata[0]['station']]
   fmidata = fulldata[0]

   # If data is fetched from multiple stations, merge the data
   if len(fulldata)>1:
      print("Merging data")
      for i in range(1,len(fulldata)):
         fmidata = mergeFMIdata(fmidata,fulldata[i])

   # Remove NaNs etc. housekeeping
   fmidata = cleanFMIdata(fmidata)

   return fmidata


#####################################################################

def dump(fmidata):
   """
   Data handling example routine: dumps fetched data to stdout.
   """
   print(u"# Station: {0} ({1}) at latitude-longitude {2}".format(fmidata['station'],fmidata['fmisid'],fmidata['latlon']))
   print("# time (UTC)",end="")
   for i in range(0,fmidata['nbvars']):
      print("\t"+fmidata['varnames'][i],end="")
   print()
   for i in range(0,fmidata['length']):
      print(datetime.datetime.utcfromtimestamp(
            fmidata['times'][i]).strftime('%Y-%m-%d %H:%M'),end="")
      for j in range(0,fmidata['nbvars']):
         print("\t"+str(round(fmidata['values'][i,j],1)),end="")
      print()
   return

#####################################################################

def write(fmidata, filename):
   """
   Data handling example routine: writes fetched data to file.
   """
   with open(filename,'w') as f:
      f.write(u"# Station: {0} ({1}) at latitude-longitude {2}\n".format(
                    fmidata['station'],
                    fmidata['fmisid'],
                    fmidata['latlon'] )
             )
      f.write("# time (UTC)")
      for i in range(0,fmidata['nbvars']):
         f.write("\t"+fmidata['varnames'][i])
      f.write("\n")
      for i in range(0,fmidata['length']):
         f.write(datetime.datetime.utcfromtimestamp(
               fmidata['times'][i]).strftime('%Y/%m/%d;%H:%M:%S:000')) # COHERENS time format
         for j in range(0,fmidata['nbvars']):
            f.write("\t"+str(fmidata['values'][i,j]))
         f.write("\n")
   return

#####################################################################


########
# MAIN
########
def main():
   if len(sys.argv)!=6:
      print("Fetch FMI weather data from nearest station.")
      print("Usage: "+sys.argv[0]+" <x> <y> <timeStart> <timeEnd> <timestep>")
      print("Example: "+sys.argv[0]+" 482900 7107700 2016-07-09T00:00:00+0000 2016-07-15T23:59:59+0000 10")
      print("x,y coordinates in ETRS-TM35FIN, times in UTC, timestep in minutes.")
      sys.exit(1)
   #
   fmidata = getFMIdata(FMIWFSURL, FMIPARAMS, FMISTATIONSFILE, MAXREQLENGTH,
                        sys.argv[3], sys.argv[4], int(sys.argv[5]), float(sys.argv[1]), float(sys.argv[2]))
   dump(fmidata)

#=== END MAIN


###################
# Default settings
###################

#---variables
UTCTIME_MIN = '2016-07-09 00:00:00+0000'
UTCTIME_MAX = '2016-08-01 23:59:59+0000'
TIMESTEP = 10 # (minutes) Min & Max time must be consistent with timestep, i.e. divisible
#POINTX = 246111.0 # Sakylan Pyhajarvi
#POINTY = 6772561.0
#POINTX = 261559.0 # Karhijarvi
#POINTY = 6835970.0 # Karhijarvi
#POINTX = 632848.0 # Puruvesi
#POINTY = 6865856.0 # Puruvesi
POINTX = 334068.0 # Lappajarvi
POINTY = 7008459.0 # Lappajarvi
#POINTX = 436446.0 # Jyvasjarvi
#POINTY = 6901246.0

#---constant settings
FMIWFSURL = 'https://opendata.fmi.fi/wfs'
FMISTATIONSFILE = 'fmi_stations_noAWOS.tsv'
MAXREQLENGTH = 7 # days
FMIPARAMS = ['temperature','windspeedms','winddirection','humidity','pressure','ri_10min','totalcloudcover']
# Units: Rain intensity (ri_10min): mm/h (10 min averaged value).

if __name__ == "__main__":
   main()

def quick(): # just return data with default values for everything (defined above)
   return getFMIdata(FMIWFSURL, FMIPARAMS, FMISTATIONSFILE, MAXREQLENGTH,
                     UTCTIME_MIN, UTCTIME_MAX, TIMESTEP, POINTX, POINTY)
