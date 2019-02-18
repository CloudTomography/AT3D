module shdom_netcdf
  ! 
  ! Dummy SHDOM netcdf interface for linking without netcdf.
  !
  implicit none
  
contains

  subroutine read_property_size_netcdf(PROPFILE, NLEG, NPX, NPY, NPZ, &
                                       NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
    character(len=*) :: propfile
    INTEGER          :: NLEG
    INTEGER          :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
    REAL             :: DELX, DELY

    print *, 'Cannot read the sizes from the netcdf property file because SHDOM was compiled without the netcdf library.'
    stop
  end subroutine read_property_size_netcdf


  subroutine read_properties_netcdf(PROPFILE, NPX, NPY, NPZ, NUMPHASE,   &
                                    MAXLEG, MAXPG, MAXPGL, DELTAM, NLEG, &
                                    PROPTYPE, ZLEVELS, MAXASYM,          &
                                    TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP) 
    character(len=*) :: propfile
    INTEGER          :: NPX, NPY, NPZ, NUMPHASE
    INTEGER          :: MAXPG, MAXPGL, MAXLEG
    LOGICAL          :: DELTAM
    integer          :: NLEG
    character        :: PROPTYPE
    real             :: zLevels(:)
    real             :: MAXASYM
    real             :: TEMPP(:), EXTINCTP(:), ALBEDOP(:), LEGENP(:)
    INTEGER        :: IPHASEP(:)

    print *, 'Cannot read the netcdf property file because SHDOM was compiled without the netcdf library.'
    stop
  end subroutine read_properties_netcdf



  subroutine output_results_netcdf_par(NX, NY, NZ, NPTS, NCELLS,                            &
                                       NSH, ML,MM,NLM, NMU,NPHI,NANG, NG,                   &
                                       PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE,  &
                                       BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,           &
                                       SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,                 &
                                       SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                                       SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, CPUTIME,  &
                                       XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID,               &
                                       FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT,      &
                                       NUMOUT, OUTTYPES, OUTPARMS, OutFileNC)
    integer :: NX, NY, NZ, NPTS, NCELLS
    integer :: NSH, ML, MM, NLM, NMU, NPHI, NANG, NG
    integer :: BCFLAG, IPFLAG, MAXITER, TOTITER, IRAD, NSHOUT, NUMOUT
    logical :: deltaM
    real    :: SOLARFLUX, SOLARMU, SOLARAZ, GNDTEMP, GNDALBEDO, SKYRAD
    real    :: SOLACC, SPLITACC, SHACC, XDOMAIN, YDOMAIN, WAVELEN, CPUTIME
    character(len=*) :: SRCTYPE, SFCTYPE, UNITS, GRIDTYPE
    character(len=*) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    real    :: XGRID(:), YGRID(:), ZGRID(:), WAVENO(:)
    real    :: fluxes(:,:,:,:), fluxdiv(:,:,:), SHTERMS(:,:,:,:), radout(:)
    character(len=1) :: OUTTYPES(:)
    real    :: OUTPARMS(:,:)
    character(len=*) :: OutFileNC                    

    print *, 'Not outputting a netcdf file because SHDOM was compiled without netcdf library.'
end subroutine output_results_netcdf_par

end module shdom_netcdf
