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
! *Usrdef_Surface_Data* User-defined surface data setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! $Date: 2015-04-10 12:30:34 +0300 (pe, 10 huhti 2015) $
!
! $Revision: 841 $
!
! Description - routines are empty but need to be declared
!
! Reference -
!
! Subroutines - usrdef_surface_absgrd, usrdef_surface_relgrd,
!               usrdef_surface_data
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

SUBROUTINE usrdef_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord,ycoord)
!************************************************************************
!
! *usrdef_surface_absgrd* Define coordinate arrays of surface grid(s)
!                         with respect to model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description - empty default routine
!             - data grid is of type 2
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
!************************************************************************
!
USE datatypes
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat) :: xcoord
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat) :: ycoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Grid file id
!*ifil*      INTEGER No. of grid file
!*n1dat*     INTEGER X-dimension of data grid
!*n2dat*     INTEGER Y-dimension of data grid
!*xcoord*    REAL    X-coordinates of data grid
!*ycoord*    REAL    Y-coordinates of data grid
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_surface_absgrd

!========================================================================

SUBROUTINE usrdef_surface_relgrd(iddesc,ifil,surfgridglb,nx,ny,nonodes)
!************************************************************************
!
! *usrdef_surface_relgrd* Define relative coordinate array of surface grid(s)
!                         with respect to model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description - empty default routine
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
!************************************************************************
!
USE datatypes
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, nonodes, nx, ny
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(nx,ny,nonodes) :: surfgridglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Grid file id
!*ifil*        INTEGER No. of grid file
!*surfgridglb* DERIVED Relative coordinates of model grid to data grid or data
!                      grid to model grid
!*nx*          INTEGER X-dimension of data grid
!*ny*          INTEGER Y-dimension of data grid 
!*nonodes*     INTEGER Number of nodes
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_surface_relgrd

!========================================================================

SUBROUTINE usrdef_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                             & novars)
!************************************************************************
!
! *usrdef_surface_data* Define surface input data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description - empty default routine
!
! Reference -
!
! Calling program - define_surface_data
!
!************************************************************************
!
USE iopars
USE syspars
USE inout_routines, ONLY: open_filepars
USE timepars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, novars, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat,novars) :: surdata

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*surdata*     REAL     Data array
!*n1dat*       INTEGER  X-dimension of data array
!*n2dat*       INTEGER  Y-dimension of data array
!*novars*      INTEGER  Number of data parameters
!
!------------------------------------------------------------------------------
!
INTEGER :: iunit
REAL :: windspeed, winddir, windn, winde
REAL :: airtemp, atmpres
REAL :: humidity, cloudcover, precipitation
REAL :: radiation
CHARACTER(LEN=2048) :: linebuffer ! Buffer for line reads

! --- First call: open file, sets %iostat to 1
IF (modfiles(iddesc,ifil,1)%iostat.EQ.0) THEN
   CALL open_filepars(modfiles(iddesc,ifil,1))
   modfiles(iddesc,ifil,1)%iostat = 1
   GOTO 1000
ENDIF

iunit = modfiles(iddesc,ifil,1)%iunit  ! Get file unit
! --- Read data
! --- Get new rows until a non-empty and non-commented row is found
readloop: DO
   READ(iunit,FMT="(A)",END=99) linebuffer
   IF ( LEN_TRIM(linebuffer).NE.0 .AND. SCAN(linebuffer,"#").NE.1 ) EXIT
ENDDO readloop

! --- Parse weather data
READ(linebuffer(1:23),FMT="(A23)") ciodatetime(1:23)
READ(linebuffer(24:),*) airtemp, windspeed, winddir, humidity, atmpres, precipitation, cloudcover

winde = 0.0
windn = 0.0
IF (windspeed>0.0) THEN
   winde = -sin(winddir*pi/180.0)*windspeed
   windn = -cos(winddir*pi/180.0)*windspeed
ENDIF

surdata(:,:,1) = winde                   ! wind speed in X (east)-direction [m/s]
surdata(:,:,2) = windn                   ! wind speed in Y (north)-direction [m/s]
!surdata(:,:,3) = atmpres*100.0          ! atmospheric pressure (Pa)
surdata(:,:,3) = airtemp                 ! air temperature [deg C]
surdata(:,:,4) = humidity/100.0          ! relative humidity 0..1
surdata(:,:,5) = MIN(cloudcover/8.0,1.0) ! cloud cover, 0..1
surdata(:,:,6) = precipitation/(60.0*60.0) ! precipitation rate [kg/m^2/s] or [mm/s]

1000 RETURN

99 modfiles(iddesc,ifil,1)%iostat = 2

END SUBROUTINE usrdef_surface_data
