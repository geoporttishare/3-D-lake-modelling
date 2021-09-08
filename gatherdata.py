#
# Gathers and processes relevant data from multiple apis and files
# and builds hydrodynamic lake model code as configured by the user.
#
# Basically just needs one coordinate within the target lake,
# model grid resolution and desired time interval.
#
# Required data sources:
# - NLS terrain database: (local) files divided into map sheets
# - Map sheets: (local) sheet shapefile (utm25LR.shp)
# - Lake shore polygons: (local) SYKE Ranta10 shapefile (jarvi10.shp)
# - Channel data: (local) SYKE channel shapefile (Uoma10.shp)
# - Traficom depth data: (remote) Open WFS service
# - Weather data: (remote) FMI's open WFS service
# - Hydrological database: (remote) SYKE's OData API

# Janne Ropponen / SYKE / oGIIR project / 2019
# Requirements and working environment:
# 1) CSC Taito cluster with geo-env module loaded
# 2) NLS terrain geodatabase
# 3) SYKE hydrological API
# 4) SYKE channel network geodatabase
# 5) Map sheet geodatabase for figuring out how NLS data is divided between files

# Last changed: 2021-09-08 / JR

from __future__ import print_function
from scipy import interpolate
from shapely.geometry import Point
from shapely.ops import unary_union

import geopandas
#import matplotlib.pyplot as plt
import numpy
import pandas
import sys
import shapely

import config # Parameters etc
import TRAFICOMdata
import SYKEdata
import FMIwfs
import ranta10
import build_model_code

###################################################################
# getBoundaryFlows
###################################################################
def getBoundaryFlows(gdf_channels, poly_lake):
   """Determines lake inflow and outflow coordinates based on where
   where channel data intersects the lake polygon borders.

   Assumes that the channel coordinates are always in outflowing direction,
   i.e. first point is an upstream point and last point is a downstream point.

   Parameters:
   channel data
   lake polygon

   Returns:
   points_outflow (list of shapely.Point)
      Coordinates of outflows
   points_inflow (list of shapely.Point) 
      Coordinates of inflows
   """
   points_inflow = []
   points_outflow = []
   # Build lists of inflow and outflow points
   for i, c in gdf_channels.iterrows():
      p0 = Point(c.geometry.xy[0][0], c.geometry.xy[1][0]) # Start coordinate of channel
      p1 = Point(c.geometry.xy[0][-1], c.geometry.xy[1][-1]) # End coordinate of channel
      #
      if (not p0.intersects(poly_lake)) and p1.intersects(poly_lake): # Inflow channel
         points_inflow.append(c.geometry.intersection(poly_lake.exterior))
      elif p0.intersects(poly_lake) and (not p1.intersects(poly_lake)): # Outflow channel
         points_outflow.append(c.geometry.intersection(poly_lake.exterior))
      elif p0.intersects(poly_lake.exterior) and p1.intersects(poly_lake.exterior): # Inflow and outflow
         points_inflow.append(p0)
         points_outflow.append(p1)
   return points_inflow, points_outflow

#-------------
# PolygonToXYZ
#-------------
def PolygonToXYZ(poly, elev=0.0):
   """Converts Polygons to list of coordinate tuples (x,y,z).

   Parameters:
   poly: shapely.Polygon
   elev: float
      Optional z level.

   Returns:
   List of polygon border coordinate tuples.
   """
   print("Converting polygon to points:")
   points = []
   # Exterior points
   x, y = poly.exterior.xy
   for p in range(1,len(x)): # Ignore first (or last) coordinate pair as it is repeated
      points.append((x[p],y[p],elev))
   print("...found",len(x),"exterior points.")
   # Interior points
   nbi = len(poly.interiors)
   print("...found", nbi, "islands.")
   for i in range(0,nbi):
      x, y = poly.interiors[i].xy
      for p in range(1,len(x)): # Ignore first (or last) coordinate pair as it is repeated
         points.append((x[p],y[p],elev))
   print("...found a total of",len(points),"shoreline points.")
   return points

###################################################################
# FeaturesToXYZ
###################################################################
def FeaturesToXYZ(geodataframe, attr, adjustment=0.0, mult=1.0):
   """Converts Linestring collections to points

   Parameters:
   GeoDataFrame geodataframe
      Object to convert
   string attr
      Attribute to convert
   float adjustment
      Optional bias value
   float mult
      Optional multiplier for unit conversions

   Returns:
   List of (x,y,z) tuples representing the Linestring coordinates.
   """
   print("Converting features to points:")
   points = []
   print("...found",len(geodataframe),"features.")
   for row in range(0,len(geodataframe)):
      x, y = geodataframe.iloc[row].geometry.xy
      z = float(geodataframe[attr].iloc[row])
      for p in range(0,max(len(x)-1,1)): # Ignore last (or first) coordinate pair as it is repeated
         points.append((x[p],y[p],mult*z+adjustment))
   print("...found a total of",len(points),"points.")
   return points

#------------------------------------------------------------
# Generates a 3-D array of depths from a collection of points
#------------------------------------------------------------
def generate_grid_griddata(points, gridresolution, interpolationmethod):
   # Generate x and y coordinates based on what area the input points cover and wanted grid resolution
   coordsx = numpy.arange(min(points[:,0])+gridresolution/2.0, max(points[:,0]), gridresolution)
   coordsy = numpy.arange(min(points[:,1])+gridresolution/2.0, max(points[:,1]), gridresolution)
   nx = len(coordsx)
   ny = len(coordsy)
   #
   # Generate a 1-D list of x, y coordinates where to interpolate point data
   coordgrid = numpy.meshgrid(coordsx,coordsy)
   coordlist = numpy.stack( (coordgrid[0].flatten(), coordgrid[1].flatten()) ).transpose()
   #
   # Simple data interpolation to grid
   interpolated = interpolate.griddata( points[:,0:2], points[:,2], coordlist, method=interpolationmethod )
   nearest = interpolate.griddata( points[:,0:2], points[:,2], coordlist, method='nearest' )
   interpolated = numpy.where(numpy.isnan(interpolated), nearest, interpolated) # Fill nans with data from nearest neighbours
   #
   # Reshape coordinate and interpolated z data into a 3-D array
   grid = numpy.stack( (coordlist[:,0],coordlist[:,1],interpolated), axis=1)
   grid3d = grid.reshape(nx, ny, 3)
   #
   return nx, ny, grid3d

#------------------------------------------------------------
# Generates a 3-D array of depths from a collection of points
#------------------------------------------------------------
def generate_grid_rbf(points, gridresolution, interpolationmethod):
   # Generate x and y coordinates based on what area the input points cover and wanted grid resolution
   coordsx = numpy.arange(min(points[:,0])+gridresolution/2.0, max(points[:,0]), gridresolution)
   coordsy = numpy.arange(min(points[:,1])+gridresolution/2.0, max(points[:,1]), gridresolution)
   nx = len(coordsx)
   ny = len(coordsy)
   xi, yi = numpy.meshgrid(coordsx,coordsy)
   #
   # Simple data interpolation to grid
   rbf = interpolate.Rbf(points[:,0],points[:,1],points[:,2])#, function=interpolationmethod, smooth=0 )
   zi = rbf(xi, yi)
   #
   # Reshape coordinate and interpolated z data into a 3-D array
   grid = numpy.stack( (xi, yi, zi), axis=1)
   grid3d = grid.reshape(nx, ny, 3)
   #
   return nx, ny, grid3d

#------------------------------------------------------------
# Generates a 3-D array of depths from a collection of points
#------------------------------------------------------------
def generate_grid_smart(points, resolution, MWL, interpolationmethod='linear'):
   """
   Smart grid interpolator that first combines all available height data points 
   to grid cells while determining the most sensible way to use the data for the
   grid cell (interpolating, averaging or selecting the most sensible value).
   Then the grid is filled by interpolating.

   The purpose is to produce a sensible calculation grid for 3-D lake modelling 
   purposes.
   """
   hix = numpy.arange(min(points[:,0]), max(points[:,0]), resolution)
   hiy = numpy.arange(min(points[:,1]), max(points[:,1]), resolution)
   nx = len(hix)
   ny = len(hiy)
   # Sort all data by y-coordinate
   sortedpoints = numpy.sort(points.view('f8,f8,f8'), order=['f1'], axis=0).view(float)
   nearest = numpy.full((nx-1,ny-1),numpy.nan)
   dbgsum = 0
   # Decimate data to regular grid cell by cell. Average original data to cells.
   # Relies heavily on the fact that the arrays being operated on are sorted correctly.
   # Assumes sortedpoints doesn't contain any points outside hix, hiy limits (especially lower)
   for j in range(0,ny-1): # row loop from bottom of the grid to the top
      #print("Row",j+1,"of",ny-1,"(",100*j/(ny-1),"% )")
      rowpoints = sortedpoints[sortedpoints[:,1]<hiy[j+1]] # Get points in this row
      rowpoints = numpy.sort(rowpoints.view('f8,f8,f8'), order=['f0'], axis=0).view(float) # Sort by columns
      sortedpoints = sortedpoints[sortedpoints[:,1]>=hiy[j+1]]
      idx = 0
      plen = len(rowpoints)
      for i in range(0,nx-1): # column loop, left to right
         localpoints = numpy.empty((0,3))
         for s in range(idx,plen):
            if rowpoints[s,0]>=hix[i+1]: # Skip to next cell if point is outside current cell
               idx = s # Start with this point at next cell
               break
            else: # Append current point to this cell
               localpoints=numpy.append(localpoints, [rowpoints[s,:]], axis=0)
               idx = s + 1 # Start with next point next time
         # Selecting cell depth value (interpolation, averaging or minimum)
         if len(localpoints)>1: # Interpolate cell value if cell contains more than 1 point
            dbgsum += len(localpoints)
            if numpy.any(localpoints[:,2]<MWL) and numpy.any(localpoints[:,2]>=MWL):
               # Border case near shorelines: select the deepest point to represent whole cell
               nearest[i,j] = min(localpoints[:,2])
            else:
               # Interpolate at center and each corner
               interpolationpoints = numpy.array( [[hix[i]+resolution/2.0, hiy[j]+resolution/2.0],
                                                   [hix[i], hiy[j]],
                                                   [hix[i], hiy[j+1]],
                                                   [hix[i+1], hiy[j]],
                                                   [hix[i+1], hiy[j+1]]
                                                  ] )
               result = interpolate.griddata(localpoints[:,0:2],localpoints[:,2],
                                             interpolationpoints,method='nearest')
               # Average points with most weight at center point
               nearest[i,j] = 0.5*result[0]+0.125*(result[1]+result[2]+result[3]+result[4])
         elif len(localpoints)==1:
            dbgsum += 1
            nearest[i,j] = localpoints[0,2]
   print("Debug: Used", dbgsum, "of", len(points), "points.")
   #
   # INTERPOLATE NaNs from regular grid
   # https://stackoverflow.com/questions/6518811/interpolate-nan-values-in-a-numpy-array
   print("Interpolating grid...")
   ix, iy = numpy.indices(nearest.shape)
   interp = numpy.array(nearest)
   interp[numpy.isnan(interp)] = interpolate.griddata(
         (ix[~numpy.isnan(nearest)], iy[~numpy.isnan(nearest)]), # points we know
         nearest[~numpy.isnan(nearest)],                    # values we know
         (ix[numpy.isnan(nearest)], iy[numpy.isnan(nearest)]),method=interpolationmethod)   # points to interpolate
   # Convert from index grids to Cartesian 3-d grid with coordinates
   hix2 = numpy.arange(min(points[:,0])+resolution/2.0, max(points[:,0])-resolution/2.0, resolution)
   hiy2 = numpy.arange(min(points[:,1])+resolution/2.0, max(points[:,1])-resolution/2.0, resolution)
   nx2 = len(hix2)
   ny2 = len(hiy2)
   xi,yi = numpy.meshgrid(hix2,hiy2,indexing='ij')
   grid3d = numpy.stack( (xi,yi,interp), axis=2 )
   #
   return nx2, ny2, grid3d

################################
# findLakeName
################################
def findLakeName(gdf, lakepoly):
   """
   Tries to determine the name of the lake polygon. Not perfect.
   Current implementation: use the name that is closest to the
   lake polygon.
   To do: make this smarter or get the name some other way.
   """
   if (len(gdf)==0):
      return ""
   distancelist = gdf.geometry.distance(lakepoly)
   mindistidx = distancelist.values.argmin()
   name = gdf.iloc[mindistidx]['TEKSTI']
   return name

###############################
# findLakeMWL
###############################
def findLakeMWL(gdf, lakepoly):
   """
   Finds lake mean water level attribute in the GeoDataFrame closest to lakepoly. 
   Not a foolproof implementation.
   """
   if (len(gdf)==0):
      return None
   mwl = -1.0
   mindist = -1.0
   for i, text in enumerate(gdf['TEKSTI'].values):
      try:
         mwltemp = float(text)
         dist = gdf.iloc[i].geometry.distance(lakepoly)
         if (mindist<0.0) or (dist<mindist):
            mindist = dist
            mwl = mwltemp
      except ValueError: # Skip erroneus mean water levels, e.g. controlled levels (127.1-128.5)
         pass
   if mwl<0.0:
      print("Fatal error, mean water level not found.")
      sys.exit(1)
   return mwl

#################################################
# cleanGrid
#################################################
def cleanGrid(grid3d, poly_lake, mwl, minheight):
   """Make the calculated grid usable in modelling:
   - remove points under mwl from outside the lake area
   - connect all wet areas and/or remove unconnected parts
   - low-pass filtering
   - minimum depth
   INCOMPLETE! TO DO: Make this better.
   """
   nx = grid3d.shape[0]
   ny = grid3d.shape[1]
   gridres = abs(grid3d[nx-1,0,0]-grid3d[0,0,0])/(nx-1) # get grid resolution from first dimension
   for i in range(0,nx):
      for j in range(0,ny):
         p = grid3d[i,j,:]
         P = Point(p[0],p[1])
         z = p[2]
         if z > minheight+50.0:
            grid3d[i,j,2] = mwl+50.0 # limit terrain height to +50 m over mwl
         elif z >= mwl and P.within(poly_lake):
            print("Warning: point", p, "within lake higher than mean water level.")
            grid3d[i,j,2] = mwl+10.0
         elif z < minheight and P.distance(poly_lake) > gridres:
            grid3d[i,j,2] = minheight # Make all land points at least minheight high
   return grid3d

###################################################
# attachFlowLocations
###################################################
def attachFlowLocations(locs, inpoints, outpoints):
   """
   Add flow location/coordinate metadata to locs dictionary.
   """
   allpoints = inpoints + outpoints
   if len(allpoints)<=0:
      print("Warning: Zero open boundaries to process.")
      return # nothing to do!
   # Assign flow point with minimum distance to each data site
   for site in locs: 
      mindist = (sys.float_info.max,-1) # (distance, index)
      for index, point in enumerate(allpoints,0):
         dist = point.distance(Point(locs[site]['x'],locs[site]['y']))
         if dist<mindist[0]:
            mindist = (dist,index)
      locs[site]['moved_x'] = allpoints[mindist[1]].x
      locs[site]['moved_y'] = allpoints[mindist[1]].y
      locs[site]['moved_distance'] = mindist[0]
      if mindist[1]<len(inpoints):
         locs[site]['type'] = 'in'
      else:
         locs[site]['type'] = 'out'
      locs[site]['enabled'] = True # Enable location by default
   return locs

##############################################
# getGridShape
##############################################
def getGridShape(grid3d, GRIDRESOLUTION, mwl):
   """Builds a list of rectangular polygons representing each
   grid cell and the exterior of the grid shape. Used in further
   processing.
   """
   nx = grid3d.shape[0]
   ny = grid3d.shape[1]
   hr = GRIDRESOLUTION/2.0
   # Build small polygon for each grid point
   polylist = []
   for i in range(0,nx):
      for j in range(0,ny):
         if ((not numpy.isnan(grid3d[i,j,2])) and grid3d[i,j,2]<mwl): # Use only wet points
            p1 = (grid3d[i,j,0]-hr,grid3d[i,j,1]-hr) # lower left
            p2 = (grid3d[i,j,0]+hr,grid3d[i,j,1]-hr) # lower right
            p3 = (grid3d[i,j,0]+hr,grid3d[i,j,1]+hr) # top right
            p4 = (grid3d[i,j,0]-hr,grid3d[i,j,1]+hr) # top left
            rect = shapely.geometry.Polygon([p1,p2,p3,p4])
            rect.i = i # add grid indices
            rect.j = j
            rect.z = grid3d[i,j,2]
            rect.c = Point((grid3d[i,j,0],grid3d[i,j,1])) # center point
            polylist.append(rect)
            #plt.plot(*polylist[-1].exterior.xy) # Debug
   #plt.show() # Debug
   # Combine polygons in polylist
   mergedpoly = shapely.ops.unary_union(polylist)
   # Try to make a single combined polygon instead of multi-part
   # TO DO: This still needs work to be usable in all cases.
   if mergedpoly.type == 'MultiPolygon':
      print("WARNING: Merging MultiPolygon forcibly by discarding orphaned areas.")
      polyareas = []
      for i, poly in enumerate(mergedpoly):
         polyareas.append( [i, poly.area] )
      polyareas.sort(key=lambda x: x[1])
      mergedpoly = mergedpoly[polyareas[-1][0]] # Keep only polygon with largest area
      #plt.plot(*mergedpoly.exterior.xy) # Debug
      #plt.show() # Debug

# Removed old merging method / JR 2019-12-12 -- might be incompatible with interface finding algorithm
#      print("MultiPolygon detected, attempting to merge by expanding.")
#      epsilon = 0.0000000001
#      while (epsilon<hr/10.0 and mergedpoly.type=='MultiPolygon'):
#         mergedpoly = mergedpoly.buffer(epsilon)
#         epsilon *= 10.0
#      if mergedpoly.type == 'MultiPolygon':
#         print("WARNING: Failed to combine MultiPolygon!")
#      else:
#         print("Success! Expanded polygons by", epsilon, " meters.")

   return polylist, mergedpoly.exterior

##########################################################
# addIntersections
##########################################################
def addIntersections(locs, polylist, gridexterior, model):
   # Finds the open boundary cell indices.
   #
   # First find the intersection coordinates of shortest distance line 
   # to grid polygon border for each flow point.
   # TODO: 1) solve ambiguities in e.g. cases where the point is equidistant to multiple
   # points on the grid border, 2) if point is inside an island (interior ring) then
   # this will still find the coordinates at the outer exterior.
   for site in locs:
      # Find coordinates of closest point to flow point at the grid exterior. For algorithm, 
      # see https://stackoverflow.com/questions/33311616/find-coordinate-of-closest-point-on-polygon-shapely
      projdist = gridexterior.project(Point(locs[site]['x'],locs[site]['y']))
      extpoint = gridexterior.interpolate(projdist)
      extcoords = extpoint.coords[0]
      px = extcoords[0]
      py = extcoords[1]
      locs[site]['grid_x'] = px # x
      locs[site]['grid_y'] = py # y
      EPSILON = 1.0E-6 # Numerical accuracy of intersection location vs polyline
      iflist = []
      for pol in polylist:
         # Check that pol shares a border with grid exterior. Excludes inner
         # polygons and polygons that touch the grid exterior only with corners.
         isect = pol.buffer(EPSILON).intersection(gridexterior)
         if ( (type(isect) is shapely.geometry.MultiLineString or
               type(isect) is shapely.geometry.LineString) and
               isect.length>EPSILON*10.0 ):
            minx = min(pol.exterior.xy[0])
            maxx = max(pol.exterior.xy[0])
            miny = min(pol.exterior.xy[1])
            maxy = max(pol.exterior.xy[1])
            hx = 0.5*(maxx-minx)+EPSILON
            hy = 0.5*(maxy-miny)+EPSILON
            if (abs(px-minx)<EPSILON and (abs(py-maxy)<hy or abs(py-miny)<hy)):
               iflist.append(['U',pol.i,pol.j,pol.z]) # U-interface at left face
            elif (abs(px-maxx)<EPSILON and (abs(py-maxy)<hy or abs(py-miny)<hy)):
               iflist.append(['U',pol.i+1,pol.j,pol.z]) # U-interface at right face
            if (abs(py-miny)<EPSILON and (abs(px-maxx)<hx or abs(px-minx)<hx)):
               iflist.append(['V',pol.i,pol.j,pol.z]) # V-interface at bottom face
            elif (abs(py-maxy)<EPSILON and (abs(px-maxx)<hx or abs(px-minx)<hx)):
               iflist.append(['V',pol.i,pol.j+1,pol.z]) # V-interface at top face
      if len(iflist)==0:
         print("Interface list is empty! This should not happen.")
         sys.exit(1)
      elif len(iflist)>1:
         # select which interface to use in case found interface matches e.g. corners
         iflist.sort(key=lambda x: x[3]) # Sort by depth (deepest first)
      locs[site]['if_indices'] = (iflist[0][1],iflist[0][2])
      locs[site]['if_orientation'] = iflist[0][0]
      if iflist[0][0] == 'U': 
         model['nrvbu'] += 1
      elif iflist[0][0] == 'V':
         model['nrvbv'] += 1
   return locs, model

#################################
# filterLocations
#################################
def filterLocations(locs, model):
   """
   Leaves only one measurement station enabled per grid location in the model
   configuration. Chooses the one with original coordinates closest 
   to the grid location.
   """
   # Possible to do list:
   # - Is there a need to check if water level data matches MWL?
   # - Should we combine multiple stations at same grid location?
   for site in locs:
      for site2 in locs:
         if (locs[site]['enabled'] and locs[site2]['enabled'] and
             locs[site]['if_indices'][0]==locs[site2]['if_indices'][0] and 
             locs[site]['if_indices'][1]==locs[site2]['if_indices'][1]):
            if locs[site]['moved_distance']<locs[site2]['moved_distance']:
               locs[site2]['enabled'] = False
               if locs[site2]['if_orientation']=='U':
                  model['nrvbu'] -= 1
               elif locs[site2]['if_orientation']=='V':
                  model['nrvbv'] -= 1
   return locs, model


#####################################
# combineEnabledLocations
#####################################
def combineEnabledLocations(W_locs, Q_locs):
   """Makes a combined dictionary of enabled W and Q locations
   in u- and v-orientations.
   """
   uvlocations = {}
   uvlocations['u'] = []
   uvlocations['v'] = []
   # Levels
   for site in W_locs:
      if W_locs[site]['enabled']:
         if W_locs[site]['if_orientation'] == 'U': # west-east if
            uvlocations['u'].append(W_locs[site])
            uvlocations['u'][-1]['id'] = site
         else: #north-south if
            uvlocations['v'].append(W_locs[site])
            uvlocations['v'][-1]['id'] = site
   # Discharges
   for site in Q_locs:
      if Q_locs[site]['enabled']:
         if Q_locs[site]['if_orientation'] == 'U': # west-east if
            uvlocations['u'].append(Q_locs[site])
            uvlocations['u'][-1]['id'] = site
         else: #north-south if
            uvlocations['v'].append(Q_locs[site])
            uvlocations['v'][-1]['id'] = site
   return uvlocations

#############################################
# getTopoDBData
#############################################
def getTopoDBData(sheets, path, fileprefix, filesuffix, classnumber, attribute):
   """Load data from all from map sheets for wanted classnumber and column.
   Returns GeoDataFrame with the requested data.
   """
   gdf = geopandas.GeoDataFrame()
   for sheet in sheets:
      print("Gathering class", classnumber, "data from map sheet", sheet)
      zipfileuri = "zip://"+path+"/"+sheet[0:2]+"/"+sheet[0:3]+"/"+sheet+".shp.zip"
      shapefilename = fileprefix+"_"+sheet+"_"+filesuffix
      gdf_all = geopandas.read_file(zipfileuri,layer=shapefilename)
      gdf = pandas.concat( [ gdf, gdf_all[gdf_all['LUOKKA']==classnumber][[attribute,'geometry']] ], axis = 0)
   return gdf

def xgetTopoDBData(sheets, path, fileprefix, filesuffix, classnumber, attribute):
   """Load data from all from map sheets for wanted classnumber and column.
   Returns GeoDataFrame with the requested data.
   """
   gdf = geopandas.GeoDataFrame()
   for sheet in sheets:
      print("Gathering class", classnumber, "data from map sheet", sheet)
      filename = path+"/"+sheet[0:2]+"/"+sheet[0:3]+"/"+fileprefix+"_"+sheet+"_"+filesuffix+".shp"
      gdf_all = geopandas.read_file(filename)
      gdf = pandas.concat( [ gdf, gdf_all[gdf_all['LUOKKA']==classnumber][[attribute,'geometry']] ], axis = 0)
   return gdf

#############################################
# writeGrid
#############################################
def writeGrid(filename, nx, ny, grid3d):
   with open(filename,'w') as f:
      for i in range(0,nx):
         for j in range(0,ny):
            p = grid3d[i,j,:]
            if numpy.isnan(p[2]):
               p[2] = -99999
            f.write(str(p[0])+"\t"+str(p[1])+"\t"+str(p[2])+"\n")
   print("Saved grid to", filename)
   return

################################################################################################
# MAIN
################################################################################################

# Get lake polygon and map sheets
poly_lake, sheets = ranta10.getLake(config)

# Load data from the map sheets (heights + mean water level + lake name)
gdf_heights     = getTopoDBData(sheets, config.TERRAINDBDIR, "k", "v", 52100, "KORARV")
gdf_depths      = getTopoDBData(sheets, config.TERRAINDBDIR, "k", "v", 54100, "KORARV")
gdf_depthpoints = getTopoDBData(sheets, config.TERRAINDBDIR, "k", "t", 54210, "TEKSTI")
gdf_MWL         = getTopoDBData(sheets, config.TERRAINDBDIR, "m", "t", 36291, "TEKSTI")
gdf_name        = getTopoDBData(sheets, config.TERRAINDBDIR, "m", "t", 36201, "TEKSTI")

# Load channel/waterway data
gdf_channels    = geopandas.read_file(config.CHANNELDB)

# Get only relevant data: depths within lake, heights surrounding lake and on islands,
# and MWL values as well as names within search radius of lake.
print("Collecting relevant data in the neighbourhood of the lake from SYKE and MML sources.")
print("...heights")
poly_larger = poly_lake.convex_hull.buffer(config.BUFFERZONE) # Area to search, slightly larger than the lake itself
gdf_heights = gdf_heights[gdf_heights.intersects(poly_larger)] # Height LineStrings

print("...depths")
gdf_depths = gdf_depths[gdf_depths.intersects(poly_lake)] # Depth LineStrings
gdf_depthpoints = gdf_depthpoints[gdf_depthpoints.intersects(poly_lake)] # Depth Points

print("...metadata")
gdf_MWL = gdf_MWL[gdf_MWL.intersects(poly_lake.buffer(config.SEARCHRADIUS))]
gdf_name = gdf_name[gdf_name.intersects(poly_lake.buffer(config.SEARCHRADIUS))]
lakename = findLakeName(gdf_name, poly_lake)
mwl = findLakeMWL(gdf_MWL, poly_lake)

print("...inflows and outflows")
gdf_channels = gdf_channels[gdf_channels.intersects(poly_lake.exterior)] # Channels that intersect lake polygon exterior
points_inflow, points_outflow = getBoundaryFlows(gdf_channels, poly_lake)

# Information
print("Detected lake", lakename, "with mean water level", mwl, "meters and area", poly_lake.area/1.0E6, "km2." )
print("Data collection area covers", len(sheets), "map sheets.")
print("The lake has", len(points_inflow), "inflows and", len(points_outflow), "outflow.")

# Make point clouds out of the depth/height data gathered
print("Generating point cloud")
depthpoints = FeaturesToXYZ(gdf_depthpoints, 'TEKSTI', mwl, -1.0)
depthpoints += FeaturesToXYZ(gdf_depths, 'KORARV', mwl, -1.0)

# If no depth data found, try to get Traficom's data for the area
if len(depthpoints) == 0:
   print("Zero depth points found in terrain database.")
   print("Checking if Finnish Transport and Communications Agency has depth data for this area.")
   depthpoints = TRAFICOMdata.getTraficomDepths(config, poly_lake, mwl)

heightpoints = FeaturesToXYZ(gdf_heights, 'KORARV')
shorepoints = PolygonToXYZ(poly_lake, mwl)

pointcloud = numpy.array( depthpoints + heightpoints + shorepoints )
pointcloud = numpy.unique(pointcloud.round(decimals=2),axis=0) # Round all coordinates to two decimals and only retain unique values

# Information
print("Point cloud has", len(depthpoints), "depth points +", 
                         len(heightpoints), "terrain points +",
                         len(shorepoints), "shoreline points.")
print(len(depthpoints)+len(heightpoints)+len(shorepoints)-len(pointcloud), 
      "duplicate points were removed.")
print("Data area (x * y):", max(pointcloud[:,0])-min(pointcloud[:,0]), "m *", 
                            max(pointcloud[:,1])-min(pointcloud[:,1]), "m.")
print("Lowest point     :", min(pointcloud[:,2]), "m.")
print("Highest point    :", max(pointcloud[:,2]), "m.")

# Make modelling grid
print("Generating grid")
#nx, ny, grid3d = generate_grid_griddata(pointcloud, config.GRIDRESOLUTION, config.INTERPOLATIONMETHOD)
#nx, ny, grid3d = generate_grid_rbf(pointcloud, config.GRIDRESOLUTION, config.INTERPOLATIONMETHOD)
nx, ny, grid3d = generate_grid_smart(pointcloud, config.GRIDRESOLUTION, mwl, config.INTERPOLATIONMETHOD)

print("Grid: nc =", nx, ", nr =", ny, ", resolution =", config.GRIDRESOLUTION, "m")
print("Grid: lowest point =", numpy.nanmin(grid3d[:,:,2]),
   "m, grid highest point =", numpy.nanmax(grid3d[:,:,2]), "m." )

# Remove points below mwl outside lake shores etc housekeeping
print("Note: skipping grid cleaning")
#print("Cleaning grid")
# Add: * Remove lake points outside lake poly
#grid3d = cleanGrid(grid3d, poly_lake, mwl, mwl+20.0)
#print("Grid: lowest point =", grid3d[:,:,2].min(), "m, grid highest point =", grid3d[:,:,2].max(), "m." )

# Calculate grid exterior polygons for later use with open boundaries
polylist, gridexterior = getGridShape(grid3d, config.GRIDRESOLUTION, mwl)

# Save grid
writeGrid(config.COHERENS_GRIDFILE, nx, ny, grid3d)

# Save point cloud
with open(config.POINTCLOUDFILE,'w') as f:
   for p in pointcloud:
      f.write(str(p[0])+"\t"+str(p[1])+"\t"+str(p[2])+"\n")
print("Saved point cloud to", config.POINTCLOUDFILE)

# Fetch and process hydrology data
print("Fetching open boundary data from SYKE hydrology API.")
W_locs, Q_locs, W_data, Q_data = SYKEdata.getHydrologicalData( config, 
                                 poly_lake.buffer(config.SEARCHRADIUS) )

# Find which inflow/outflow is attached to each hydrology data. Add keys: moved_x, moved_y, type, enabled
W_locs = attachFlowLocations(W_locs, points_inflow, points_outflow)
Q_locs = attachFlowLocations(Q_locs, points_inflow, points_outflow)

# Find grid indices where each outflow/inflow should be located and in which direction U/V. 
# Add keys: grid_i, grid_j
model = {}
model['nrvbu'] = 0
model['nrvbv'] = 0
W_locs, model = addIntersections(W_locs, polylist, gridexterior, model)
Q_locs, model = addIntersections(Q_locs, polylist, gridexterior, model)

# Disable unusable locations (enabled=False) and update model configuration.
W_locs, model = filterLocations(W_locs, model)
Q_locs, model = filterLocations(Q_locs, model)

# Make a combined dictionary of enabled locations
uvlocations = combineEnabledLocations(W_locs, Q_locs)

# Save open boundary data
SYKEdata.writeOBData(uvlocations, W_data, Q_data)

# Fetch and process weather data from FMI open WFS
fmidata = FMIwfs.getFMIdata(config.FMIWFSURL, config.FMIPARAMS, 
                            config.FMISTATIONSFILE, config.FMIMAXREQLENGTH,
                            config.UTCTIME_MIN, config.UTCTIME_MAX, 
                            config.FMITIMESTEP, config.POINTX, config.POINTY)

FMIwfs.write(fmidata,config.COHERENS_WEATHEROUTFILE)

# Build COHERENS model code
model = build_model_code.build(uvlocations, W_locs, model, nx, ny, mwl, grid3d)

print("All done!")
# TO DO: make flags to turn features on or off, e.g. whether to gather time series data
