#
# COHERENS user definitions model code builder
# oGIIR project
#
# Janne Ropponen, Finnish Environment Institute
# 2019-2021
#
# Changelog:
# ----------
# 2021-09-08 Pyproj changes
# 2019-09-20 First version

from string import Template
from dateutil import parser as dateparser
from fractions import Fraction

import codecs
import pyproj

import config

def coherenstimestring(s):
   """
   Returns time string s in COHERENS model compatible format: YYYY/MM/DD;hh:mm:ss:fff
   """
   ns = s[0:4]   + "/" + s[5:7]   +"/" + s[8:10]  + ";" \
      + s[11:13] + ":" + s[14:16] +":" + s[17:19] + ":" + s[20:23]
   return ns

def build(locslist, W_locs, model, nx, ny, mwl, grid3d):
   """
   Build COHERENS hydrodynamic model compatible user definition code.
   Generates open boundary indices based on locations and their types,
   grid definitions and timestep information and writes basic 
   Usrdef_Model.f90 code.
   """
   # Generate COHERENS code for model indices
   print("Generating COHERENS model code")

   modtypes = {}
   modtypes['transport'] = "4" # ityp2dobu/ityp2dobv
   modtypes['elevation'] = "3" # ityp2dobu/ityp2dobv
   modtypes['elevationfile'] = "2" # iobc2dtype: file contains elevation data
   modtypes['dischargefile'] = "3" # iobc2dtype: file contains depth-integrated current data
   
   model['iobu'] = ""
   model['jobu'] = ""
   model['iobv'] = ""
   model['jobv'] = ""
   model['ityp2dobu'] = ""
   model['ityp2dobv'] = ""
   model['iobc2dtype'] = ""
   
   # Build Usrdef_Model code for Open boundaries (U)
   for l in range(0,model['nrvbu']):
      # ob to grid index mapping (usrdef_grid)
      model['iobu'] += "iobu(" + str(l+1) + ") = " \
                       + str(locslist['u'][l]['if_indices'][0]+1) + " ! " \
                       + locslist['u'][l]['name'] + " (" \
                       + locslist['u'][l]['type'] + ")\n" # iobu(l) = i
      model['jobu'] += "jobu(" + str(l+1) + ") = " \
                       + str(locslist['u'][l]['if_indices'][1]+1) + " ! " \
                       + locslist['u'][l]['name'] + " (" \
                       + locslist['u'][l]['type'] + ")\n" # jobu(l) = j
      # ob type mapping (usrdef_2dobc_spec)
      if locslist['u'][l]['id'] in W_locs: # current ob found in elevation list
         model['ityp2dobu'] += "ityp2dobu(" + str(l+1) + ") = " + modtypes['elevation'] + "\n" # ityp2dobu(l) = type
         model['iobc2dtype'] += "iobc2dtype(" + str(l+2) + ") = " + modtypes['elevationfile'] + "\n" # iobc2dtype
      else: # current ob in discharge list (Q_locs)
         model['ityp2dobu'] += "ityp2dobu(" + str(l+1) + ") = " + modtypes['transport'] + "\n" # ityp2dobu(l) = type
         model['iobc2dtype'] += "iobc2dtype(" + str(l+2) + ") = " + modtypes['dischargefile'] + "\n" # iobc2dtype
   
   # Open boundaries (V)
   for l in range(model['nrvbu'],model['nrvbu']+model['nrvbv']):
      # ob to grid index mapping (usrdef_grid)
      model['iobv'] += "iobv(" + str(l+1) + ") = " \
                       + str(locslist['v'][l-model['nrvbu']]['if_indices'][0]+1) + " ! " \
                       + locslist['v'][l-model['nrvbu']]['name'] + " (" \
                       + locslist['v'][l-model['nrvbu']]['type'] + ")\n" # iobv(l) = ii
      model['jobv'] += "jobv(" + str(l+1) + ") = " \
                       + str(locslist['v'][l-model['nrvbu']]['if_indices'][1]+1) + " ! " \
                       + locslist['v'][l-model['nrvbu']]['name'] + " (" \
                       + locslist['v'][l-model['nrvbu']]['type'] + ")\n" # iobv(l) = jj
      # ob type mapping (usrdef_2dobc_spec)
      if locslist['v'][l-model['nrvbu']]['id'] in W_locs: # current ob found in elevation list
         model['ityp2dobv'] += "ityp2dobv(" + str(l-model['nrvbu']+1) + ") = " + modtypes['elevation'] + "\n" # ityp2dobv(l) = type
         model['iobc2dtype'] += "iobc2dtype(" + str(l+2) + ") = " + modtypes['elevationfile'] + "\n" # iobc2dtype
      else: # current ob in discharge list (Q_locs)
         model['ityp2dobv'] += "ityp2dobv(" + str(l-model['nrvbu']+1) + ") = " + modtypes['transport'] + "\n" # ityp2dobv(l) = type
         model['iobc2dtype'] += "iobc2dtype(" + str(l+2) + ") = " + modtypes['dischargefile'] + "\n" # iobc2dtype
   
   # Model flags
   model['CStartDateTime'] = coherenstimestring(config.UTCTIME_MIN[0:23])
   model['CEndDateTime'] = coherenstimestring(config.UTCTIME_MAX[0:23])
   dt0 = dateparser.parse(config.UTCTIME_MIN) # convert to datetime object
   dt1 = dateparser.parse(config.UTCTIME_MAX)
   runlength = (dt1-dt0).total_seconds() # Run length in seconds

   # TO DO: Smarter 2d timestep selection
   delt2d = config.GRIDRESOLUTION/50.0
   if delt2d>config.COHERENS_MAXDELT2D:
      delt2d = config.COHERENS_MAXDELT2D
   elif config.COHERENS_MAXDELT2D%delt2d>config.COHERENS_EPSILON:
      delt2d = config.COHERENS_MAXDELT2D/math.ceil(config.COHERENS_MAXDELT2D/delt2d)

   model['delt2d'] = delt2d
   
   # TO DO: Smart 3d timestep selection
   model['ic3d'] = config.COHERENS_IC3D
   
   #time_zone = 2.0
   model['time_zone'] = 0.0 # UTC
   
   transformer = pyproj.Transformer.from_crs('epsg:3067','epsg:4258') # etrstm35fin -> etsr89 (~wgs84)
   coords = transformer.transform(config.POINTX,config.POINTY)
   model['dlat_ref'] = coords[0]
   model['dlon_ref'] = coords[1]

   model['nc'] = nx
   model['nr'] = ny
   model['nz'] = config.COHERENS_NZ
   
   model['delxdat'] = config.GRIDRESOLUTION
   model['delydat'] = config.GRIDRESOLUTION
   model['x0dat'] = grid3d[0,0,0]-config.GRIDRESOLUTION/2.0
   model['y0dat'] = grid3d[0,0,1]-config.GRIDRESOLUTION/2.0
   
   model['meanlevel'] = mwl
   
   model['grdfile'] = config.COHERENS_GRIDFILE
   model['metfile'] = config.COHERENS_WEATHEROUTFILE
   
   # Substitute data from template
   with codecs.open(config.COHERENS_USRDEFMODEL_TEMPLATE,mode='r',encoding='utf-8') as file:
      template = Template(file.read())
   
   usrdef_model_code = template.safe_substitute(model)
   
   # Write code to Usrdef_Model.f90
   with codecs.open(config.COHERENS_USRDEFMODEL,mode='w',encoding='utf-8') as file:
      file.write(usrdef_model_code)

   return model
