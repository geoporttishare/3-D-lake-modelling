!
! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!
! http://ec.europa.eu/idabc/eupl
!
! Unless required by the applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.

!************************************************************************
!
! *Usrdef_Model* User-defined model setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! $Date: 2015-04-10 12:30:34 +0300 (pe, 10 huhti 2015) $
!
! $Revision: 841 $
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Routines - usrdef_init_params, usrdef_mod_params, usrdef_grid,
!            usrdef_partition, usrdef_phsics, usrdef_1dsur_spec,
!            usrdef_2dobc_spec, usrdef_profobc_spec, usrdef_1dsur_data,
!            usrdef_2dobc_data, usrdef_profobc_data, usrdef_rlxobc_spec
!
!************************************************************************
!

!*******************************************************
! Modified for use in automatic lake model building code
!
! oGIIR-project
!
! Last changed: 2019-12-13
! Janne Ropponen, Finnish Environment Institute
!*******************************************************

!========================================================================

SUBROUTINE usrdef_init_params
!************************************************************************
!
! *usrdef_init_params* Define parameters for monitoring
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - simulation_start
!
!************************************************************************
!
USE iopars
USE paralpars

!---program leveling in log files
levprocs_ini = 7
levprocs_run = 3
IF (npworld.GT.1) THEN
   levprocs_ini(2:npworld) = 0
   levprocs_run(2:npworld) = 0
ENDIF

RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params* Define model parameters
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE physpars
USE switches
USE syspars
USE timepars

INTEGER :: i
CHARACTER(LEN=100) :: obfilename

!============
! MODEL SETUP
!============

!---Global options
IF (npworld.GE.1) nprocscoh = npworld ! Detect available threads at runtime

!---Model options
iopt_adv_2D = 1 ! 2D advection scheme: 1 = Upwind (default)
iopt_adv_3D = 3 ! 3D advection scheme: 3 = TVD

iopt_dens = 2 ! Density formulation
iopt_fld = 1 ! Enable drying/wetting (1=without dynamic masks, 2=using dynamic masks)
iopt_temp = 2 ! Temperature calculation mode

iopt_meteo = 1 ! Enable surface forcing (weather)
iopt_meteo_data = 1 ! Meteo input data: uwindatc, vwindatc, atmpres, airtemp, relhum, cloud_cover
iopt_meteo_pres = 0 ! Disable atmospheric pressure input
iopt_meteo_stres = 1 ! Enable input of wind surface stress
iopt_meteo_heat = 1 ! Enable temperature/heat flux input
iopt_meteo_precip = 1 ! precipitation type: 0 = disable, 1 = precipitation, 2 = evap-precip
IF (iopt_meteo_precip .GT. 0) iopt_sflux_precip = 3 ! Salinity flux

!iopt_sed = 1 ! enable sediment/suspended solids calculation module (w/o flocculation)
!iopt_obc_sed = 1 ! enable sediment open boundaries

!---Physical parameters / freshwater lake specific settings at northern latitudes
optattcoef1_cst = 10.0 ! Infrared (10.0)
optattcoef2_cst = 2.0/8.0 ! shortwave (0.067). optattcoef = 2.0/Secchi depth
opt_frac = 0.54 ! Infrared fraction of irradiance absorbed at sea surface (0.54)
!albedo_os = 0.10 ! In Surface_Fluxes.F90 ! Surface albedo (0.06), Korhonen/Suomen vesistojen lampotilaolot 1900-luvulla: 5-10%
specheat = 4181.4 ! Specific heat of water (3987.5 J/kg/deg C)
temp_ref = 4.0
temp_min = float_fill ! Allow water temperatures below 0.0 degC
sal_ref = 0.0
zrough_cst = 0.10 ! [m] Bottom roughness

!---Time settings
CStartDateTime(1:23) = "${CStartDateTime}" ! YYYY.MM.DD.HH.MM.SS.000
CEndDateTime(1:23) = "${CEndDateTime}"
delt2d = ${delt2d} ! 2-D timestep
ic3d = ${ic3d} ! Calculate 3-D values every ic3d 2-D timestep
time_zone = ${time_zone}

!---Grid setup
dlat_ref = ${dlat_ref} ! Latitude
dlon_ref = ${dlon_ref} ! Longitude
nc = ${nc}+1 ! Grid dimension: columns (x)
nr = ${nr}+1 ! Grid dimension: rows (y)
nz = ${nz} ! Grid dimension: depth levels
!iopt_grid_vtype = 1 ! sigma coordinates
iopt_grid_vtype = 3 ! user defined vertical coordinates
iopt_grid_vtype_transf = 0 ! user defined vertical coordinates
surfacegrids(igrd_model,1)%delxdat = ${delxdat} ! Grid size x [m]
surfacegrids(igrd_model,1)%delydat = ${delydat} ! Grid size y [m]
surfacegrids(igrd_model,1)%x0dat = ${x0dat} ! Lower left corner x
surfacegrids(igrd_model,1)%y0dat = ${y0dat} ! Lower left corner y
surfacegrids(igrd_meteo,1)%nhtype = 4

!---Open boundaries
iopt_obc_2D = 1 ! enable 2D open boundaries
nrvbu = ${nrvbu} ! Number of u-direction (x) river boundaries
nrvbv = ${nrvbv} ! Number of v-direction (y) river boundaries

!--- Output options
iopt_out_tsers = 1 ! enable time series output
iopt_CDF_tlim = 1 ! NetCDF output time mode
nosetstsr = 1 ! Number of output time series file sets
novarstsr = 6 ! Number of output time series variables

!=========================
! Input/output files setup
!=========================
!modfiles(key id,file number,in/out)
! File status: %status = '0' (not activated), 'R' (read), 'W' (write), 'N' (user defined)
! File format: %form = 'A' (ASCII), 'U' (binary), 'N' (NetCDF)
! Data update: %tlims = (/ start, end, step /) (in multiples of delt2d)

!---Model grid file
modfiles(io_modgrd,1,1)%status = 'N'
modfiles(io_modgrd,1,1)%form = 'A'
modfiles(io_modgrd,1,1)%filename = "${grdfile}"

!---Open boundary condition files (rivers)
IF ((${nrvbu}+${nrvbv}).GT.0) THEN
   modfiles(io_2uvobc,1,1)%status = 'N'
   DO i=1,(${nrvbu}+${nrvbv})
      modfiles(io_2uvobc,i+1,1)%status = 'N'
      modfiles(io_2uvobc,i+1,1)%form = 'A'
      modfiles(io_2uvobc,i+1,1)%tlims = (/ 0, int_fill, ic3d /)
      WRITE(obfilename, "(A6,I0.3,A4)") "obdata", i, ".dat"
      modfiles(io_2uvobc,i+1,1)%filename = TRIM(obfilename)
   ENDDO
ENDIF

!---Surface forcing files (weather)
modfiles(io_metsur,1,1)%status = 'N'
modfiles(io_metsur,1,1)%form = 'A'
modfiles(io_metsur,1,1)%tlims = (/ 0, int_fill, ic3d /)
modfiles(io_metsur,1,1)%filename = "${metfile}"

!---Sediment specifiers and open boundary conditions files
!modfiles(io_sedspc,1,1)%status = 'N'
!modfiles(io_sedobc,1,1)%status = 'N' ! sediment obc specifiers 
!modfiles(io_sedobc,2,1)%status = 'N' ! sediment obc data

!---Initial conditions setups
modfiles(io_inicon,ics_phys,1)%status = 'N' ! Set up physical initial conditions in usrdef_phsics

RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define arrays for grid file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE iopars
USE depths
USE grid
USE gridpars
USE inout_routines, ONLY: open_filepars
USE switches

INTEGER :: i,j,k,iunit
REAL :: xcoord, ycoord, zlevel, mindep, meanlevel, mindepth2, layerheight

mindep = 0.3 ! Minimum cell depth at mean water level height
meanlevel = ${meanlevel} ! Mean water level

CALL open_filepars(modfiles(io_modgrd,1,1))
iunit=modfiles(io_modgrd,1,1)%iunit

rowloop: DO i=1,nc-1
colloop: DO j=1,nr-1
   READ(iunit,*) xcoord, ycoord, depmeanglb(i,j)
   IF (depmeanglb(i,j).LT.0.0) THEN
      depmeanglb(i,j) = meanlevel ! convert negative points to land
   ENDIF
   IF ((depmeanglb(i,j).LT.meanlevel) .AND. (depmeanglb(i,j).GT.(meanlevel-mindep))) THEN
      depmeanglb(i,j) = meanlevel-mindep ! set minimum depth to mindep
   ENDIF
   depmeanglb(i,j) = MAX(meanlevel-depmeanglb(i,j),0.0) ! reclassify heights to depths

!--- s-coordinate constant 2x top & 2x bottom layers when depth is over minimum
   layerheight = 0.50 ! top layer constant thickness
   mindepth2 = nz*layerheight ! minimum total depth for constant layer
   IF (iopt_grid_vtype.EQ.3 .AND. iopt_grid_vtype_transf.EQ.0) THEN
      IF (depmeanglb(i,j).LT.mindepth2 .AND. depmeanglb(i,j).GT.0.0) THEN
         gscoordglb(i,j,:) = (/ (REAL(k-1)*1.0/REAL(nz),k=1,nz+1) /)
      ELSEIF (depmeanglb(i,j).GT.0.0) THEN
         ! top layers
         gscoordglb(i,j,nz+1) = 1.0
         gscoordglb(i,j,nz) = (depmeanglb(i,j)-layerheight)/depmeanglb(i,j)
         gscoordglb(i,j,nz-1) = (depmeanglb(i,j)-2.0*layerheight)/depmeanglb(i,j)

         ! bottom layers
         gscoordglb(i,j,3) = 2.0*layerheight/depmeanglb(i,j)
         gscoordglb(i,j,2) = layerheight/depmeanglb(i,j)
         gscoordglb(i,j,1) = 0.0

         ! Middle layers
         gscoordglb(i,j,4:nz-2) = &
            & (/ ( gscoordglb(i,j,3)+REAL(k-3)*(gscoordglb(i,j,nz-1)-gscoordglb(i,j,3))/REAL(nz-4),k=4,nz-2) /)
 
     ENDIF
   ENDIF

ENDDO colloop
ENDDO rowloop

!---Inflow and outflow grid locations (open boundaries)
${iobu}
${jobu}
${iobv}
${jobv} 

RETURN

END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decomposition for parallel mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - domain_decomposition
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_partition

!========================================================================

SUBROUTINE usrdef_phsics
!************************************************************************
!
! *usrdef_phsics* Define arrays for initial conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE depths

! --- Set initial water level to mean water level
! TO DO: set level to actual level at starting time (if available)
zeta(:,:) = 0.0 ! ${meanlevel}

RETURN

END SUBROUTINE usrdef_phsics

!========================================================================

SUBROUTINE usrdef_1dsur_spec
!************************************************************************
!
! *usrdef_1dsur_spec* Define specifier arrays for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_1dobc_spec
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!                     at U- and V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
USE gridpars
USE iopars
USE obconds

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name          Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!
!------------------------------------------------------------------------------
!
INTEGER :: l, nofiles

nofiles = maxdatafiles(io_2uvobc,1) ! maximum number of data files

${ityp2dobu}
${ityp2dobv}
${iobc2dtype}

!*ityp2dobu*    INTEGER Type of 2-D o.b.c. at U-nodes
! (1:nobu)        =  2  => zero gradient for U, elevation specified at internal node
!                 =  3  => specified elevation and U from local solution
!                 =  4  => specified transport
!                 = 14  => specified depth mean current
!                 = 15  => specified discharge (at local position)
!                 = 16  => discharge specified along o.b. sections

!*ityp2dobv*    INTEGER Type of 2-D o.b.c. at V-nodes

!*iobc2dtype*   INTEGER Type of data in 2-D data files
! (2:nofiles)     =  1  => transports and elevations
!                 =  2  => elevations
!                 =  3  => transports

!ityp2dobu(1:)     = 3  ! River boundary: specified elevation
!ityp2dobu(:)     = 4  ! River boundary: specified transport
!iobc2dtype(2:nofiles) = 2  ! All files use elevations only

! Default open boundary to file indexing (each file maps to one open boundary)
no2dobuv(2:nofiles)    = 1 ! Every file represents data for 1 open boundary
DO l=1,nobu+nobv
   index2dobuv(l,l+1) = l ! Mapping of ob numbers to file numbers
END DO

RETURN

END SUBROUTINE usrdef_2dobc_spec

!========================================================================

SUBROUTINE usrdef_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                             & iprofrlx,noprofsd,indexprof,indexvar,novars,&
                             & nofiles,nobux,nobvy)
!************************************************************************
!
! *usrdef_profobc_spec* Define specifier arrays for open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_profobc_spec
!
!************************************************************************
!
USE gridpars
USE iopars
USE relaxation

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, nofiles, novars
INTEGER, INTENT(INOUT), DIMENSION(2:nofiles) :: noprofsd
INTEGER, INTENT(INOUT), DIMENSION(nobux) :: itypobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy) :: itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexprof
INTEGER, INTENT(INOUT), DIMENSION(novars*(nobux+nobvy),2:nofiles) :: indexvar
INTEGER, INTENT(INOUT), DIMENSION(norlxzones) :: iprofrlx

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*iprofrlx*  INTEGER Disables/enables relaxation at open boundary zones
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof* INTEGER Mapping array of the profile numbers in the data files to
!                    the profile numbers assigned to the open boundaries. The
!                    physical size of the first dimension equals the number of
!                    profiles in a data file.
!*indexvar*  INTEGER Defines the variable number of the profiles in a data
!                    file. The physical size of the first dimension equals the
!                    number of profiles in a data file.
!*novars*    INTEGER Total number of variables
!*nofiles*   INTEGER Number of data files (+1)
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
!
!------------------------------------------------------------------------------
!

!---TO DO: Implement profiles and/or sediment/tracer variables

!SELECT CASE (iddesc)
!   CASE (io_sedobc)
!      noprofsd = 1
!      iprofobvy(1,1) = 1
!      indexprof(1,2) = 1
!END SELECT

RETURN

END SUBROUTINE usrdef_profobc_spec

!========================================================================

SUBROUTINE usrdef_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *usrdef_1dsur_data* Define data for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_1dsur_data
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: novars
REAL, INTENT(INOUT), DIMENSION(novars) :: data1d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time in data file
!*data1d*      REAL    Surface data
!*novars*      INTEGER Number of surface data
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_1dsur_data

!========================================================================

SUBROUTINE usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *usrdef_2dobc_data* Define open boundary data for 2-D mode at U- and V-nodes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_2dobc_data
!
!************************************************************************
!
USE gridpars
USE iopars
USE obconds
USE syspars
USE inout_routines, ONLY: open_filepars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: data2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Data file id (io_2uvobc or io_2xyobc)
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*data2d*      REAL    Input data
!*nodat*       INTEGER Number of input data
!*novars*      INTEGER Number of input parameters
!
!------------------------------------------------------------------------------
!
INTEGER :: iunit
REAL :: discharge, gridwidth, elevation

CHARACTER(LEN=512) :: linebuffer ! Buffer for line reads

! -------- On first call open file and set status --------
IF (modfiles(io_2uvobc,ifil,1)%iostat.EQ.0 ) THEN
  CALL open_filepars(modfiles(io_2uvobc,ifil,1))   ! open file at first call
  modfiles(io_2uvobc,ifil,1)%iostat=1              ! file ready to read
  RETURN
ENDIF

! --- Select local grid cell width based on whether this file
!     corresponds to U-open boundary or V-open boundary
IF (index2dobuv(ifil-1,ifil).LE.nrvbu) THEN
   gridwidth = surfacegrids(igrd_model,1)%delydat ! U-open boundary, regular grid assumed
ELSE
   gridwidth = surfacegrids(igrd_model,1)%delxdat ! V-open boundary, regular grid assumed
ENDIF

iunit = modfiles(io_2uvobc,ifil,1)%iunit

! --- Get new row until we get a non-empty and non-commented row
readloop: DO
   READ(iunit,FMT="(A)",END=99) linebuffer
   IF ( LEN_TRIM(linebuffer).NE.0 .AND. SCAN(linebuffer,"#").NE.1 ) EXIT
ENDDO readloop

IF (iobc2dtype(ifil).EQ.3) THEN
   ! -------- Read river discharge data from file(s) --------
   READ(linebuffer(1:23),FMT="(A23)") ciodatetime(1:23)   ! Read time
   READ(linebuffer(24:),*) discharge
   data2d = discharge/gridwidth ! Convert discharge to specified transport
ELSEIF (iobc2dtype(ifil).EQ.2) THEN
   ! -------- Read river elevation data from file(s) --------
   READ(linebuffer(1:23),FMT="(A23)") ciodatetime(1:23)   ! Read time
   READ(linebuffer(24:),*) elevation
   data2d = elevation-${meanlevel} ! Lake mean level
ELSE
   ! TO DO: Add support for other o.b.c. types
   WRITE(*,*) "Warning:, iobc2dtype(", ifil, "=", iobc2dtype(ifil)
ENDIF

RETURN

99 modfiles(io_2uvobc,ifil,1)%iostat = 2

END SUBROUTINE usrdef_2dobc_data

!========================================================================

SUBROUTINE usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,numprofs)
!************************************************************************
!
! *usrdef_profobc_data* Define physical open boundary profiles
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.2
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_profobc_data
!
!************************************************************************
!
USE gridpars
USE iopars
USE syspars
USE timepars

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, numprofs
REAL, INTENT(INOUT), DIMENSION(numprofs,nz) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*psiprofdat*  REAL     Profile arrays
!*numprofs*    INTEGER  Number of profiles in data file
!
!------------------------------------------------------------------------------
!

!---TO DO: Implement profile and tracer data reading

!IF (modfiles(iddesc,2,1)%iostat.EQ.0) THEN
!   modfiles(iddesc,2,1)%iostat = 1
!   GOTO 1000
!ENDIF
!
!ciodatetime = CStartDateTime
!
!SELECT CASE (iddesc)
!CASE (io_sedobc)
!   psiprofdat(1,:) = 0.5 ! constant concentration of tracer
!END SELECT
!
!
!1000 RETURN

END SUBROUTINE usrdef_profobc_data

!========================================================================

SUBROUTINE usrdef_rlxobc_spec
!************************************************************************
!
! *usrdef_rlxobc_spec* Define specifier arrays for relaxation
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_rlxobc_spec
!
!************************************************************************
!


RETURN

END SUBROUTINE usrdef_rlxobc_spec
