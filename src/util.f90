subroutine average_subpixel_rays (npixels,nrays, weighted_stokes, nstokes, &
                              pixel_index, observables)

  integer npixels, nrays, nstokes, pixel_index(nrays)
!f2py intent(in) pixel_index, npixels, nrays
  real weighted_stokes(nstokes, nrays)
!f2py intent(in) stokes
  real observables(nstokes, npixels)
!f2py intent(out) observables

  integer i,j
  double precision temp(nstokes)

  iray=1
  do i=1,npixels
    temp = 0.0D0
    do while (pixel_index(iray) + 1 .EQ. i)
      temp(:) = temp(:) + weighted_stokes(:,iray)
      iray = iray + 1
    enddo
    observables(:,i) = temp(:)
  enddo
  return
end
