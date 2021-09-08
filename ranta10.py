# Janne Ropponen / SYKE / oGIIR project / 2019
# Requirements and working environment:
# 1) CSC Taito cluster with geo-env module loaded
# 2) SYKE Ranta10 lake shapefiles

from __future__ import print_function
from shapely.geometry import Point
import copy
import geopandas
import sys
import locale # For UTF8 output stuff
import codecs # For UTF8 output stuff


# TO DO list
#-----------
# Handle multi-part lakes, e.g. Paijanne, Kemijarvi
# Algorithm?:
# 0) poly_lake = gdf_lake.iloc[0].geometry
# 1) Check if first found polygon intersects with 1 or more other
#    polygons in Ranta10 jarvet: 
#    gdf_new=gdf_lakes[gdf_lakes.intersects(poly_lake)]
# 2) Unary union of new polygons
#    poly_lake = shapely.ops.unary_union([gdf_new.iloc[i].geometry for i in len(gdf_new)])
# 3) Check again

def getLakeShape(SYKERANTA10SHP,X,Y):
   """
   Returns the lake shape that contains coordinates X, Y.
   
   Arguments:
   string     SYKERANTA10SHP Full path to lake shape data
   float      X              X coordinate inside lake in ETRS-TM35FIN
   float      Y              Y coordinate inside lake in ETRS-TM35FIN
   #
   Returns:
   shapely Polygon       Lake shape that that matches the given coordinates or None
   """
   print("Searching for lake geometry in", SYKERANTA10SHP)
   gdf_lakes = geopandas.read_file(SYKERANTA10SHP)
   P = Point(X, Y)
   #
   gdf_lake = gdf_lakes[gdf_lakes.contains(P)]
   if gdf_lake.empty or len(gdf_lake)>1:
      print("Ambiguous point.", len(gdf_lake), "matches found.")
      sys.exit(1)  
   print("Extracted lake", gdf_lake.iloc[0]['nimi'], "with area", 
         gdf_lake.iloc[0]['Shape_STAr']/10**6, "km2.")
   return gdf_lake


def getSheets(MAPSHEETFILE, gdf, BUFFERZONE):
   """
   Extracts names of map sheets that intersect with the shape gdf.

   Arguments:
   string       MAPSHEETFILE  Full path to map sheet data shapefile
   GeoDataFrame gdf           Single polygon (len(gdf=1))
   float        BUFFERZONE    Amount of buffer zone (m) to add to gdf geometry
   
   Returns:
   list of strings            Names of map sheets that intersect with gdf
   """
   print("Selecting map sheets from",MAPSHEETFILE)
   gdf_mapsheets = geopandas.read_file(MAPSHEETFILE)

   # Select only sheets that intersect buffered lake polygon
   gdf_buffered = copy.deepcopy(gdf)
   gdf_buffered['geometry'] = gdf_buffered.convex_hull.buffer(BUFFERZONE) # Add buffer zone
   gdf_sheets = geopandas.overlay(gdf_mapsheets, 
                                  gdf, 
                                  how='intersection')
   
   if gdf_sheets.empty:
      print("Sheet list empty, aborting.")
      sys.exit(1)

   # Generate list of sheet names
   sheets = []
   for i in range(1,len(gdf_sheets)):
      sheets.append(str(gdf_sheets.iloc[i]['LEHTITUNNU']))
   print("Found",len(gdf_sheets),"map sheets:",sheets)

   return sheets


def getLake(config):
   """Determines the lake polygon and the map sheets that contain data
   for the lake based on given coordinates.
   """
   gdf_lake = getLakeShape(config.SYKERANTA10SHP,config.POINTX,config.POINTY)
   sheets = getSheets(config.MAPSHEETFILE, gdf_lake, config.BUFFERZONE)
   poly_lake = gdf_lake.iloc[0].geometry
   return poly_lake, sheets


if __name__ == "__main__":
   GRIDRESOLUTION = 500.0 # Grid resolution in meters
   # Bufferzone outside the lake within which height data is collected.
   # Also acts as a guarantee that we get correct MWL and lake name points 
   # within the loaded map sheets.
   BUFFERZONE = 100.0 + 3.0*GRIDRESOLUTION
   #
   #POINTX = 261559.0 # Karhijarvi
   #POINTY = 6835970.0 # Karhijarvi
   POINTX = 436446.0 # Jyvasjarvi
   POINTY = 6901246.0
   SYKERANTA10SHP = "/appl/data/geo/syke/ranta10jarvet/jarvi10.shp"
   MAPSHEETFILE = "/users/jroppone/ogiir/karttaruudut/utm25LR.shp"
   #
   gdf_lake = getLakeShape(SYKERANTA10SHP,POINTX,POINTY)
   #
   sheetstoload = getSheets(MAPSHEETFILE, gdf_lake, BUFFERZONE) # Get sheet names
   #
   poly_lake = gdf_lake.iloc[0].geometry
   poly_bounds = poly_lake.convex_hull.buffer(BUFFERZONE)

