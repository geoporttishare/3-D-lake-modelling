#!/usr/bin/python
#encoding: utf-8

from __future__ import print_function
from shapely.geometry import Point
from shapely.geometry import box
import sys
import json
import geopandas
import time
from owslib.wfs import WebFeatureService

##########################################################################
# Downloads Traficom's lake or sea depth data Web Feature Service (WFS).
#
# Available data: https://julkinen.traficom.fi/oskari/?zoomLevel=2&coord=463617.19864277664_6943478.691860745&mapLayers=3+100+,30+100+syvyyskayra,31+100+Syvyyspiste_point&uuid=55bc3b05-5e13-451a-bca0-a25560000fd1&noSavedState=true
#
# oGIIR project
# Janne Ropponen/SYKE
# Last changed: 2021-09-08
##########################################################################

####################################################################
# getWFSDepthPoints
####################################################################
def getWFSDepthPoints(wfsurl, xmin, ymin, xmax, ymax):
   """Gets depth data (points and curves) from the Traficom WFS

   Parameters
   ----------
   wfsurl: string
             Traficom WFS url
   xmin, ymin, xmax, ymax: float
             Bounding rectangle corner coordinates in ETRS-TM35FIN

   Returns a list of tuples.
   The tuples are coordinates of data points (x, y, depth).
   """
   ### Get depth data from Traficom WFS
   # Note, might also want to save some metadata from the wfs query, 
   # eg. measurement date and data quality. See:
   # https://vayla.fi/documents/20473/38174/Meri_tietosis%C3%A4ll%C3%B6n_kuvaus.pdf/78afa9e5-8f7c-4430-b798-f9848c79123f
   # Instead of json outputFormat, another options is to use zipped 
   # shapefile (outputFormat="shape-zip").
   wfs = WebFeatureService(url=wfsurl, version='2.0.0')
   points_json = wfs.getfeature(typename='rajoitettu:Sounding_P', 
                 bbox=(xmin,ymin,xmax,ymax), outputFormat="json") 
   lines_json = wfs.getfeature(typename='rajoitettu:DepthContour_L', 
                bbox=(xmin,ymin,xmax,ymax), outputFormat="json")
   
   ### Convert json response to dict
   try:
      d_points = json.loads(points_json.getvalue())
      d_lines = json.loads(lines_json.getvalue())
   except AttributeError:
      # Sometimes repeated identical queries fail and give a ResponseWrapper result instead of cStringIO
      print("Unable to fetch data... maybe wait a few minutes before making this exact query again.")
      return []
   
   nb_points = d_points['totalFeatures']
   nb_lines = d_lines['totalFeatures']
   print("Traficom: Found", nb_points, "depth points.")
   print("Traficom: Found", nb_lines, "depth contours.")
   
   ### Convert Point and LineString dicts to tuples of depth points (x,y,z)
   ### and combine the data
   points_tuple = []
   for i in range(0, nb_points):
      p = d_points['features'][i]['geometry']['coordinates']
      points_tuple.append( (p[0], p[1], p[2]) )
   
   lines_tuple = []
   for i in range(0, nb_lines):
      depth = d_lines['features'][i]['properties']['VALDCO'] # Depth of this depth curve
      for mls in d_lines['features'][i]['geometry']['coordinates']: # One or more LineStrings in a curve
         if isinstance(mls[0],float): # Is this the already coordinates?
            lines_tuple.append( (mls[0], mls[1], float(depth)) )
         else: # No, iterate over coordinates in LineString.
            for coords in mls: # One or more coordinates in a LineString
               lines_tuple.append( (coords[0], coords[1], float(depth)) )
   
   print("Traficom: Combining", len(lines_tuple), "line points and", 
                                len(points_tuple), "depth points.")
   
   return lines_tuple + points_tuple

########################################################
# getTraficomDepths
########################################################
def getTraficomDepths(config, poly_lake, mwl=0.0):
   """Fetch and process depth points from Traficom's WFS.
   Points are processed into a list of coordinates.
   """
   xmin, ymin, xmax, ymax = poly_lake.bounds
   depthpoints = []
   depthtuple = getWFSDepthPoints(config.TRAFICOMURL, xmin, ymin, xmax, ymax)
   if len(depthtuple) == 0:
      print("Sorry, no depth points found for this lake")
   else:
      if ( len(depthtuple)>config.MAXCOMPLEXITY and 
           len(poly_lake.exterior.coords)>config.MAXCOMPLEXITY ):
         print("Simplifying lake geometry to save calculation time")
         simplepoly = poly_lake.buffer(0.5*config.SIMPLIFYTOLERANCE).simplify(config.SIMPLIFYTOLERANCE) # To maintain lake area, first grow by half of simplifying tolerance
      else:
         simplepoly = poly_lake
      print("Selecting relevant points (this might take some time).")
      liviGS = geopandas.GeoSeries([Point(p) for p in depthtuple])
      l = len(liviGS)
      t0 = time.time()
      t = time.time()
      for i,p in enumerate(liviGS):
         if (time.time()-t)>10.0:
            t=time.time()
            print("  ", format(100.0*i/l,'.1f'), "%, ETA", format((t-t0)*(l-i)/i,'.2f'), "s")
         if p.intersects(simplepoly):
            depthpoints.append((p.x,p.y,mwl-p.z))
   return depthpoints


####################################################################
# MAIN
####################################################################
def main(config):
   if len(sys.argv)!=5:
      print("Fetches lake/marine depth point cloud from Traficom WFS service.")
      print("Usage: python "+sys.argv[0]+" <xmin> <ymin> <xmax> <ymax>")
      print("Example (Oulujarvi region): python "+sys.argv[0]+" 482900 7107700 554000 7160000")
      print("xmin etc are bounding box coordinates in ETRS-TM35FIN.")
      print("See data available at: https://www.traficom.fi/en/statistics-and-publications/spatial-dataset-material")
      sys.exit(1)

   poly_lake = box(float(sys.argv[1]), 
                   float(sys.argv[2]),
                   float(sys.argv[3]), 
                   float(sys.argv[4]))

   pointlist = getTraficomDepths(config, poly_lake)

   ### Example usage: dump points to screen.
   print("#License: https://www.traficom.fi/fi/tilastot-ja-julkaisut/paikkatietoaineistot/aineistojen-kaytto-ja-lisenssit")
   for point in pointlist:
      print(str(point[0])+"\t"+str(point[1])+"\t"+str(point[2]))
# END main

class config:
   TRAFICOMURL = "https://julkinen.traficom.fi/inspirepalvelu/rajoitettu/wfs?request=getcapabilities"
   MAXCOMPLEXITY = 1000 # Maximum lake complexity before simplifying geometry (for performance purposes)
   SIMPLIFYTOLERANCE = 100.0 # [m] Simplifying tolerance.

if __name__ == "__main__":
   main(config)

