! Fortran subroutines written by
! Jesse Loveridge, University of Illinois at Urbana-Champaign, 2020-2021
! for inclusion in pyshdom. `average_subpixel_rays` is essential to the operation
! of pyshdom.
! `cell_average` is useful for diagnostics when comparing between grids.
! There is also a ray integration routine to allow linear tomography methods
! utilizing an SHDOM-type grid.


subroutine wigner_transform(p1, p2, p3, p4, coef, nrank, wts, nangle, &
  mu)
! This subroutine converts a 4 component phase matrix as a function of angle
! to the wigner function coefficient expansion.
! This is useful for converting from tabulations as a function of angle
! to phase functions that can be used in SHDOM.
! It does not include the possible p5, p6 entries from non-spherical particles.
  implicit none

  integer nangle, nrank
!f2py intent(in) :: nangle, nrank
  double precision wts(nangle), coef(6,0:nrank), mu(nangle)
!f2py intent(in) :: wts, mu
!f2py intent(out) :: coef
  double precision p1(nangle), p2(nangle), p3(nangle), p4(nangle)
!f2py intent(in) :: p1, p2, p3, p4

  double precision D00(0:nrank), D20(0:nrank),D22P(0:nrank)
  double precision D22M(0:nrank)
  double precision a2, a3, f
  integer i,l

  coef = 0.0d0

  do i=1, nangle
    call wignerfct(mu(i), nrank, 0,0, D00)
    call wignerfct(mu(i), nrank, 2,2, D22P)
    call wignerfct(-mu(i), nrank, 2,2, D22M)
    call wignerfct(mu(i), nrank, 2,0, D20)

    do l=0, nrank
      f = (-1.0d0)**l
      coef(1,l) = coef(1,l) + (l+0.5D0)*wts(i)*P1(i)*D00(l)
      coef(2,l) = coef(2,l) + (l+0.5D0)*wts(i)*(P1(i)+P3(i))*D22P(l)
      coef(3,l) = coef(3,l) + (l+0.5D0)*wts(i)*(P1(i)-P3(i))*f*D22M(l)
      coef(4,l) = coef(4,l) + (l+0.5D0)*wts(i)*P3(i)*D00(l)
      coef(5,l) = coef(5,l) - (l+0.5D0)*wts(i)*P2(i)*D20(l)
      coef(6,l) = coef(6,l) - (l+0.5D0)*wts(i)*P4(i)*D20(l)
    enddo
  enddo

  do l = 0, nrank
    a2 = 0.5d0*(coef(2,l) + coef(3,l))
    a3 = 0.5d0*(coef(2,l) - coef(3,l))
    coef(2,l) = a2
    coef(3,l) = a3
  enddo

end subroutine






subroutine grid_smoothing(xgrid, ygrid, zgrid, nx, ny, nz, &
  field, cost, gradient, weights, mode, ierr, errmsg, &
  direction_weights, huber_parameter)
! This subroutine evaluates volume integrals of different norms of 3D
! derivatives of `field` assuming field is defined on grid points
! with coordinates given by `xgrid`, `ygrid`, `zgrid` and a linear
! interpolation kernel is assumed for the norm of the gradient.
! Norms on the gradient field may be L2 or Pseudo-Huber, which is a differentiable
! form of the L1 norm. An 'L1' norm may be specified but this is not differentiable
! around zero which can lead to artifacts.
! This subroutine is useful for supplying regularizing constraints on fields
! in an inverse problem.
  implicit none
  integer nx, ny, nz
  double precision xgrid(nx), ygrid(ny), zgrid(nz), field(nz,ny,nx)
  double precision weights(nz,ny,nx)
  double precision gradient(nz,ny,nx), cost
  double precision direction_weights(3)
  double precision huber_parameter
  character(len=2) mode
  integer ierr
  character(len=600) errmsg
!f2py intent(in) :: nx, ny, nz, xgrid, ygrid, zgrid, field, weights
!f2py intent(in) :: direction_weights, huber_parameter
!f2py intent(out) :: gradient, cost, ierr, errmsg
  integer i,j,k
  double precision dx,dy,dz,invdx,invdy,invdz,a,b,c
  double precision a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4
  double precision x_derivs, y_derivs, z_derivs, volume
  double precision norm_of_increment, cost_increment

  cost = 0.0D0
  gradient = 0.0D0
  ierr = 0
  a = direction_weights(1)
  b = direction_weights(2)
  c = direction_weights(3)

  ! for each cell calculate the volume integral of
  ! the square of each derivative.
  do i = 1, nx-1
    dx = xgrid(i+1) - xgrid(i)
    invdx = 1.0D0/dx

    do j = 1, ny-1
      dy = ygrid(j+1) - ygrid(j)
      invdy = 1.0D0/dy

      do k = 1, nz-1
        dz = zgrid(k+1) - zgrid(k)
        invdz = 1.0D0/dz

        ! each directional derivative can vary linearly across
        ! the remaining plane.
        volume = dx*dy*dz*0.125D0*(weights(k,j,i) + weights(k,j,i+1) &
                 + weights(k,j+1,i) + weights(k,j+1,i+1) &
                 + weights(k+1,j,i+1)+weights(k+1,j+1,i) &
                 + weights(k+1,j+1,i+1) + weights(k+1,j,i))
        x_derivs = 0.25D0*invdx*((field(k,j,i+1) - field(k,j,i))&
                 + (field(k,j+1,i+1) - field(k,j+1,i)) &
                 + (field(k+1,j,i+1) - field(k+1,j,i)) &
                 + (field(k+1,j+1,i+1) - field(k+1,j+1,i)))
        y_derivs = 0.25D0*invdy*((field(k,j+1,i) - field(k,j,i))&
                 + (field(k,j+1,i+1) - field(k,j,i+1)) &
                 + (field(k+1,j+1,i) - field(k+1,j,i)) &
                 + (field(k+1,j+1,i+1) - field(k+1,j,i+1)))
        z_derivs = 0.25D0*invdz*((field(k+1,j,i) - field(k,j,i))&
                 + (field(k+1,j,i+1) - field(k,j,i+1)) &
                 + (field(k+1,j+1,i) - field(k,j+1,i)) &
                 + (field(k+1,j+1,i+1) - field(k,j+1,i+1)))

        if (mode .eq. 'l1') then
          cost = cost + volume*(a*abs(x_derivs) + b*abs(y_derivs) + &
                                c*abs(z_derivs))
!         This gives the sign function for minimization of the l1 norm.
!         Note that this is non-smooth near zero.
!         If this is bad then I will add a smoothed region with a hyperparameter
!         which will reduce to quadratic behavior for small gradients.
!         Note this also assumes that direction_weights are positive.

          a1 = sign(a,field(k,j,i+1) - field(k,j,i))
          a2 = sign(a,field(k,j+1,i+1) - field(k,j+1,i))
          a3 = sign(a,field(k+1,j,i+1) - field(k+1,j,i))
          a4 = sign(a,field(k+1,j+1,i+1) - field(k+1,j+1,i))

          b1 = sign(b,field(k,j+1,i) - field(k,j,i))
          b2 = sign(b,field(k,j+1,i+1) - field(k,j,i+1))
          b3 = sign(b,field(k+1,j+1,i) - field(k+1,j,i))
          b4 = sign(b,field(k+1,j+1,i+1) - field(k+1,j,i+1))

          c1 = sign(c,field(k+1,j,i) - field(k,j,i))
          c2 = sign(c,field(k+1,j,i+1) - field(k,j,i+1))
          c3 = sign(c,field(k+1,j+1,i) - field(k,j+1,i))
          c4 = sign(c,field(k+1,j+1,i+1) - field(k,j+1,i+1))

          gradient(k,j,i) = gradient(k,j,i) + volume*(    &
            -a1*invdx - b1*invdy - c1*invdz)
          gradient(k,j,i+1) = gradient(k,j,i+1) + volume*(    &
            a1*invdx - b2*invdy - c2*invdz)
          gradient(k,j+1,i) = gradient(k,j+1,i) + volume*(    &
            -a2*invdx +b1*invdy - c3*invdz)
          gradient(k,j+1,i+1) = gradient(k,j+1,i+1) + volume*(    &
            a2*invdx +b2*invdy - c4*invdz)
          gradient(k+1,j,i) = gradient(k+1,j,i) + volume*(    &
            -a3*invdx - b3*invdy + c1*invdz)
          gradient(k+1,j,i+1) = gradient(k+1,j,i+1) + volume*(    &
            a3*invdx - b4*invdy + c2*invdz)
          gradient(k+1,j+1,i) = gradient(k+1,j+1,i) + volume*(    &
            -a4*invdx +b3*invdy + c3*invdz)
          gradient(k+1,j+1,i+1) = gradient(k+1,j+1,i+1) + volume*(    &
            a4*invdx +b4*invdy + c4*invdz)

!          print *, k,j,i, gradient(:,1,1)

        else if (mode .eq. 'l2') then
          cost = cost + volume*(a*x_derivs*x_derivs + b*y_derivs*y_derivs + &
                                c*z_derivs*z_derivs)

          gradient(k,j,i) = gradient(k,j,i) + volume*(    &
            -a*x_derivs*invdx - b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k,j,i+1) = gradient(k,j,i+1) + volume*(    &
            a*x_derivs*invdx - b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k,j+1,i) = gradient(k,j+1,i) + volume*(    &
            -a*x_derivs*invdx +b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k,j+1,i+1) = gradient(k,j+1,i+1) + volume*(    &
            a*x_derivs*invdx +b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k+1,j,i) = gradient(k+1,j,i) + volume*(    &
            -a*x_derivs*invdx - b*y_derivs*invdy + c*z_derivs*invdz)
          gradient(k+1,j,i+1) = gradient(k+1,j,i+1) + volume*(    &
            a*x_derivs*invdx - b*y_derivs*invdy + c*z_derivs*invdz)
          gradient(k+1,j+1,i) = gradient(k+1,j+1,i) + volume*(    &
            -a*x_derivs*invdx +b*y_derivs*invdy + c*z_derivs*invdz)
          gradient(k+1,j+1,i+1) = gradient(k+1,j+1,i+1) + volume*(    &
            a*x_derivs*invdx +b*y_derivs*invdy + c*z_derivs*invdz)
        else if (mode .eq. 'ph') then
          norm_of_increment = a*x_derivs*x_derivs + b*y_derivs*y_derivs + c*z_derivs*z_derivs
          cost_increment = huber_parameter**2 * (sqrt(1 + norm_of_increment/(huber_parameter**2)) - 1)

          ! Update the volume variable with the huber weight to the derivative.
          cost = cost + cost_increment*volume
          volume = volume/sqrt(1 + norm_of_increment/(huber_parameter**2))
          ! The derivative is then just the l2 derivative with the Huber
          ! weight added.
          gradient(k,j,i) = gradient(k,j,i) + volume*(    &
            -a*x_derivs*invdx - b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k,j,i+1) = gradient(k,j,i+1) + volume*(    &
            a*x_derivs*invdx - b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k,j+1,i) = gradient(k,j+1,i) + volume*(    &
            -a*x_derivs*invdx +b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k,j+1,i+1) = gradient(k,j+1,i+1) + volume*(    &
            a*x_derivs*invdx +b*y_derivs*invdy - c*z_derivs*invdz)
          gradient(k+1,j,i) = gradient(k+1,j,i) + volume*(    &
            -a*x_derivs*invdx - b*y_derivs*invdy + c*z_derivs*invdz)
          gradient(k+1,j,i+1) = gradient(k+1,j,i+1) + volume*(    &
            a*x_derivs*invdx - b*y_derivs*invdy + c*z_derivs*invdz)
          gradient(k+1,j+1,i) = gradient(k+1,j+1,i) + volume*(    &
            -a*x_derivs*invdx +b*y_derivs*invdy + c*z_derivs*invdz)
          gradient(k+1,j+1,i+1) = gradient(k+1,j+1,i+1) + volume*(    &
            a*x_derivs*invdx +b*y_derivs*invdy + c*z_derivs*invdz)
        else
          ierr = 1
          write(errmsg, *) "GRID_SMOOTHING: Unrecognized argument &
                            for 'mode'", mode
        endif
      enddo
    enddo
  enddo

  gradient = 0.25D0*gradient
  if (mode .eq. 'l2') then
    gradient = gradient*2
  endif

end subroutine grid_smoothing

subroutine phase_function_mixing(scatter_coefficients, phase_tables, &
  phase_indices, interpolation_weights, &
  nparticles, maxpg, maxleg, nstleg, mixed_phase_table, mixed_phase_indices, &
  asym_tol, phase_tol, dolp_tol, numphase, maxnphase, nangles, maxnmicro, &
  nphase, ierr, errmsg)
! This subroutine combines optical properties using similar logic to SHDOM's
! PROPGEN program (see propgen.f90). Adapted for pyshdom by Jesse Loveridge.
! A greedy algorithm is used to select the set of phase functions:
!   New phase functions are added to the set if either the
!   maximum relative error across nangles of the phase function or maximum relative
!   error in the asymmetry parameter are exceeded across all of the set of
!   phase functions.
!
!   We set a hard upper limit on the total number of phase functions
!   to use and the method will fail if exceeded. Full arrays are allocated
!   based on largest possible sizes rather than writing to file.

  implicit none
  integer maxpg, maxleg, nstleg, nparticles, numphase, maxnmicro
  real scatter_coefficients(maxpg, nparticles)
!f2py intent(in) :: scatter_coefficients
  real interpolation_weights(maxnmicro, maxpg, nparticles)
!f2py intent(in) :: interpolation_weights
  real phase_tables(nstleg, 0:maxleg, numphase)
!f2py intent(in) :: phase_tables
  integer phase_indices(maxnmicro,maxpg, nparticles)
!f2py intent(in) :: phase_indices
  double precision asym_tol, phase_tol, dolp_tol
!f2py intent(in) :: asym_tol, phase_tol, dolp_tol
  integer maxnphase, nangles
!f2py intent(in) :: maxphase, nangles
  real mixed_phase_table(nstleg, 0:maxleg, maxnphase)
!f2py intent(out) :: mixed_phase_table
  integer mixed_phase_indices(maxpg)
!f2py intent(out) :: mixed_phase_indices
  integer ierr, nphase
  character(len=600) errmsg
!f2py intent(out) :: ierr, errmsg, nphase

  integer i, l, jpart, k, minind
  real asym(maxnphase), phase(2,nangles, maxnphase)
  real asym0, phase0(2,nangles), minerr, err
  real dolp(nangles), dolp0(nangles)
  double precision scatter
  real wigcoef(nstleg, 0:maxleg)

  ierr = 0

  mixed_phase_table = 0.0
  asym = 0.0
  phase = 0.0
  mixed_phase_indices = 0
  nphase = 0

  do i=1,maxpg
    wigcoef = 0.0
    scatter = 0.0
    do jpart=1,nparticles
      scatter = scatter + scatter_coefficients(i,jpart)
      do k=1,maxnmicro
        wigcoef = wigcoef + scatter_coefficients(i,jpart)* &
          interpolation_weights(k,i,jpart)* &
          phase_tables(:,:,phase_indices(k,i,jpart))
      enddo
    enddo
    if (scatter > 0.0) then

      wigcoef = wigcoef / scatter
      asym0 = wigcoef(1, 1)/3
      ! this is in shdomsub5.f
      call phasefunc_from_wigner(nangles, maxleg, wigcoef(:,0:), phase0, &
        2, ierr, errmsg, .True., 6)

      if (nphase .eq. 0) then
      ! initialize with the first phase function we find.
        nphase = nphase + 1
        asym(nphase) = asym0
        phase(:,:,nphase) = phase0
        mixed_phase_table(:,:,nphase) = wigcoef
        mixed_phase_indices(i) = nphase
      else
      ! compare with existing phase functions
        minind = 1
        minerr = 1000.0
        do l=1,nphase
          dolp0 = phase0(2,:)/phase0(1,:)
          dolp = phase(2,:,l)/phase(1,:,l)
          err = abs(asym0-asym(l))/asym_tol &
            + maxval(abs((phase0(1,:) - phase(1,:,l))/max(0.001,phase(1,:,l))))/ &
            phase_tol + maxval(abs((dolp0(:) - dolp(:))/ &
              max(0.001, dolp(:))))/dolp_tol

          if (err < minerr) then
            minerr = err
            minind = l
          endif
        enddo
        ! if smallest error is sufficiently small then use that phase function
        ! otherwise add this phase function to the set of phase functions.
        dolp = phase(2,:,minind)/phase(1,:,minind)
        if (abs(asym0 - asym(minind)) < asym_tol .and. phase_tol > &
          maxval(abs((phase0(1,:)-phase(1,:,minind))/max(0.001,phase(1,:,minind)))) &
          .and. dolp_tol > &
          maxval(abs((dolp0(:)-dolp)/max(0.001,dolp))) ) then
          mixed_phase_indices(i) = minind
        else
          nphase = nphase + 1
          if (nphase > maxnphase) then
            ierr = 1
            write(errmsg, *) 'Maximum number of phase functions exceeded when ', &
              'mixing optical properties of different particle species. Either ', &
              'decrease accuracy tolerance or increase maxnphase.', &
              ' nphase=',nphase, ' maxnphase=', maxnphase, ' gridpoint=', i, &
              ' total_number_of_gridpoints=', maxpg
            return
          endif
          mixed_phase_table(:,:,nphase) = wigcoef
          asym(nphase) = asym0
          phase(:,:,nphase) = phase0
          mixed_phase_indices(i) = nphase
        endif
      endif
    endif
  enddo
  return
end subroutine phase_function_mixing

subroutine adjoint_linear_interpolation(xgrid, ygrid, zgrid, &
  nx,ny,nz, output, inxgrid, inygrid, inzgrid, innx,inny, &
  innz, field)
! This subroutine is written for mapping between two grids with
! linear interpolation in 3D.
! If the unknowns are on a coarser resolution
! grid and linearly interpolated to the property grid for use in SHDOM
! we need to correctly perform the adjoint of the interpolation
! to transform the gradients from the property grid back up to the
! unknowns. This subroutine performs this adjoint interpolation.
! This subroutine assumes that input has been sanitized and does
! no error checking.

  implicit none
  integer nx,ny,nz,innx,inny, innz
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
! f2py intent(in) :: xgrid, ygrid, zgrid
  double precision inxgrid(innx), inygrid(inny), inzgrid(innz)
! f2py intent(in) :: inxgrid, inygrid, inzgrid
  double precision output(nx,ny,nz), field(innx,inny,innz)
!f2py intent(out) :: output
!f2py intent(in) :: field

  integer i,j,k
  do k=1,innz
    do j=1,inny
      do i=1,innx
        call adjoint_interp_point(xgrid, ygrid, zgrid, &
         nx,ny,nz, field(i,j,k), inxgrid(i), inygrid(j), &
         inzgrid(k), output)
      enddo
    enddo
  enddo

  return
end subroutine adjoint_linear_interpolation


subroutine adjoint_interp_point(xgrid, ygrid, &
  zgrid, nx,ny,nz, fieldval, x,y, z, output)

  implicit none
  integer nx, ny, nz
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
  double precision x,y,z, fieldval
  double precision output(nx,ny,nz)

  double precision u,v,w,f1,f2,f3,f4,f5,f6
  double precision f7,f8
  integer il, iu, im, ix, iy, iz
! binary searches for the position.
  il=0
  iu=nx
  do while (iu-il .gt. 1)
    im = (iu+il)/2
    if (x .ge. xgrid(im)) then
      il = im
    else
      iu=im
    endif
  enddo
  ix = max(il,1)
  u = (x - xgrid(ix))/(xgrid(ix+1) - xgrid(ix))
  u = max(min(u, 1.0D0), 0.0D0)

  il=0
  iu=ny
  do while (iu-il .gt. 1)
    im = (iu+il)/2
    if (x .ge. ygrid(im)) then
      il = im
    else
      iu=im
    endif
  enddo
  iy = max(il,1)
  v = (y - ygrid(iy))/(ygrid(iy+1) - ygrid(iy))
  v = max(min(v, 1.0D0), 0.0D0)

  il=0
  iu=nz
  do while (iu-il .gt. 1)
    im = (iu+il)/2
    if (x .ge. zgrid(im)) then
      il = im
    else
      iu=im
    endif
  enddo
  iz = max(il,1)
  w = (z - zgrid(iz))/(zgrid(iz+1) - zgrid(iz))
  w = max(min(w, 1.0D0), 0.0D0)
! calculate interpolation weights.
  f1 = (1-u)*(1-v)*(1-w)
  f2 =    u *(1-v)*(1-w)
  f3 = (1-u)*   v *(1-w)
  f4 =    u *   v *(1-w)
  f5 = (1-u)*(1-v)*   w
  f6 =    u *(1-v)*   w
  f7 = (1-u)*   v *   w
  f8 =    u *   v *   w
! add to the output field the contributions of the
! adjoint interpolation.
  output(ix,iy,iz) = output(ix,iy,iz) + f1*fieldval
  output(ix+1,iy,iz) = output(ix+1,iy,iz) + f2*fieldval
  output(ix,iy+1,iz) = output(ix,iy+1,iz) + f3*fieldval
  output(ix+1,iy+1,iz) = output(ix+1,iy+1,iz) + f4*fieldval
  output(ix,iy,iz+1) = output(ix,iy,iz+1) + f5*fieldval
  output(ix+1,iy,iz+1) = output(ix+1,iy,iz+1) + f6*fieldval
  output(ix,iy+1,iz+1) = output(ix,iy+1,iz+1) + f7*fieldval
  output(ix+1,iy+1,iz+1) = output(ix+1,iy+1,iz+1) + f8*fieldval

  return
end subroutine adjoint_interp_point



subroutine average_subpixel_rays (npixels,nrays, weighted_stokes, nstokes, &
                              pixel_index, observables)
! Averages over sub-pixel rays to calculate pixel average observables.
! See pyshdom.containers._calculate_observables.

  implicit none
  integer npixels, nrays, nstokes, pixel_index(nrays)
!f2py intent(in) pixel_index, npixels, nrays
  real weighted_stokes(nstokes, nrays)
!f2py intent(in) stokes
  real observables(nstokes, npixels)
!f2py intent(out) observables

  integer i,j, iray, pixind
  double precision temp(nstokes)

  iray=1
  pixind = 0
  do i=1,npixels
    temp = 0.0D0
    do while (pixind + 1 .eq. i)
      temp(:) = temp(:) + weighted_stokes(:,iray)
      iray = iray + 1
      if (iray .ge. nrays) then
        ! this breaks the loop at iray=nray-1
        pixind = i+100
      else
        pixind = pixel_index(iray)
      endif
    enddo
    observables(:,i) = temp(:)
  enddo
! add last ray to the last pixel. Very lazy way to do things.
  observables(:,npixels) = observables(:,npixels) + weighted_stokes(:,nrays)
  return
end subroutine average_subpixel_rays

subroutine nearest_binary(size, reference, val, index_low, &
  index_high)

  implicit none
  integer size, index_low, index_high, upper, lower
!f2py intent(out) index_low, index_high
  double precision reference(size), val
!f2py intent(in) reference, val
  logical go
  integer index_test

  upper = size
  lower = 1
  go = .true.

  do while (go)
    index_test = (upper + lower) / 2 !rounds down.
    if (reference(index_test) .eq. val) then
      index_low = index_test
      index_high = index_test + 1
      go = .false.
    else if (reference(index_test) > val) then
      upper = index_test
    elseif (reference(index_test) < val) then
      lower = index_test
    endif

    if (upper - lower .eq. 1) then
      index_low = lower
      index_high = upper
      go = .false.
    endif
    enddo

  return
end subroutine nearest_binary

subroutine cell_average(nx1,ny1,nz1, &
  nx2,ny2,nz2, values1, values2, xgrid1,ygrid1,zgrid1, &
  xgrid2,ygrid2,zgrid2, ref_vol, ref_val, other_vol, other_val)

  double precision xgrid1(nx1),ygrid1(ny1),zgrid1(nz1)
!f2py intent(in) xgrid1, ygrid1, zgrid1
  double precision xgrid2(nx2),ygrid2(ny2),zgrid2(nz2)
!f2py intent(in) xgrid2, ygrid2, zgrid2
  double precision values1(nx1,ny1,nz1), values2(nx2,ny2,nz2)
!f2py intent(in) values1, values2
  double precision ref_vol(nx1-1,ny1-1,nz1-1), ref_val(nx1-1,ny1-1,nz1-1)
!f2py intent(out) ref_vol, ref_val
  double precision other_vol(nx1-1,ny1-1,nz1-1), other_val(nx1-1,ny1-1,nz1-1)
!f2py intent(out) other_vol, other_val

  double precision dx1,dy1,dz1,dx2,dy2,dz2
  double precision min_x,min_y,min_z,max_x,max_y,max_z
  integer ic,jc,kc, i,j,k
  integer index_low(3), index_high(3), temp_low, temp_high
  integer gridptr(3,8), t
  double precision f(8), location(3)

  other = 0.0D0
  reference = 0.0D0
  do kc=1,nz1 - 1
    do jc=1,ny1 - 1
      do ic=1,nx1 - 1

        dx1 = xgrid1(ic+1) - xgrid1(ic)
        dy1 = ygrid1(jc+1) - ygrid1(jc)
        dz1 = zgrid1(kc+1) - zgrid1(kc)
        location(1) = xgrid1(ic) + dx1/2
        location(2) = ygrid1(jc) + dy1/2
        location(3) = zgrid1(kc) + dz1/2

        ! call util_locate_point(nx1,ny1,nz1,xgrid1,ygrid1,zgrid1, &
        !     location, gridptr)
        call util_get_interp_kernel2(ic,jc,kc, location(1), location(2), &
          location(3), f,nx1,ny1,nz1,xgrid1,ygrid1,zgrid1)

        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(1)* &
                                  values1(ic,jc,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(2)* &
                                  values1(ic+1,jc,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(3)* &
                                  values1(ic,jc+1,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(4)* &
                                  values1(ic+1,jc+1,kc)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(5)* &
                                  values1(ic,jc,kc+1)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(6)* &
                                  values1(ic+1,jc,kc+1)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(7)* &
                                  values1(ic,jc+1,kc+1)
        ref_val(ic,jc,kc) = ref_val(ic,jc,kc) +  dx1*dy1*dz1*f(8)* &
                                  values1(ic+1,jc+1,kc+1)

        ref_vol(ic,jc,kc) = dx1*dy1*dz1
        !for further speed up, replace these searches with an increment scheme
        !as index_low_new = index_high_old when ic is incremented (i think.)
        call nearest_binary(nx2, xgrid2, xgrid1(ic), index_low(1), &
          temp_high)
        call nearest_binary(nx2, xgrid2, xgrid1(ic+1), temp_low, &
          index_high(1))
        call nearest_binary(ny2, ygrid2, ygrid1(jc), index_low(2), &
          temp_high)
        call nearest_binary(ny2, ygrid2, ygrid1(jc+1), temp_low, &
          index_high(2))
        call nearest_binary(nz2, zgrid2, zgrid1(kc), index_low(3), &
          temp_high)
        call nearest_binary(nz2, zgrid2, zgrid1(kc+1), temp_low, &
          index_high(3))

        do k=index_low(3),index_high(3) - 1
            do j=index_low(2),index_high(2) - 1
              do i=index_low(1),index_high(1) - 1

                min_x = min(max(xgrid2(i), xgrid1(ic)),xgrid1(ic+1))
                min_y = min(max(ygrid2(j), ygrid1(jc)),ygrid1(jc+1))
                min_z = min(max(zgrid2(k), zgrid1(kc)),zgrid1(kc+1))

                max_x = max(min(xgrid2(i+1), xgrid1(ic+1)),xgrid1(ic))
                max_y = max(min(ygrid2(j+1), ygrid1(jc+1)),ygrid1(jc))
                max_z = max(min(zgrid2(k+1), zgrid1(kc+1)),zgrid1(kc))
                dx2 = max_x - min_x
                dy2 = max_y - min_y
                dz2 = max_z - min_z

                location(1) = 0.5*(min_x + max_x)
                location(2) = 0.5*(min_y + max_y)
                location(3) = 0.5*(min_z + max_z)

                ! call util_locate_point(nx2,ny2,nz2,xgrid2,ygrid2,zgrid2, &
                !     location, gridptr)
                call util_get_interp_kernel2(i,j,k, location(1), location(2), &
                  location(3), f,nx2,ny2,nz2,xgrid2,ygrid2,zgrid2)

                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(1)* &
                                          values2(i,j,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(2)* &
                                          values2(i+1,j,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(3)* &
                                          values2(i,j+1,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(4)* &
                                          values2(i+1,j+1,k)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(5)* &
                                          values2(i,j,k+1)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(6)* &
                                          values2(i+1,j,k+1)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(7)* &
                                          values2(i,j+1,k+1)
                other_val(ic,jc,kc) = other_val(ic,jc,kc) +  dx2*dy2*dz2*f(8)* &
                                          values2(i+1,j+1,k+1)

                other_vol(ic,jc,kc)  = other_vol(ic,jc,kc) + dx2*dy2*dz2

              enddo
            enddo
          enddo


      enddo
    enddo
  enddo
end subroutine cell_average


subroutine util_integrate_rays(nrays, nx, ny, nz, xgrid, ygrid, zgrid, &
  start_positions, end_positions, weights, paths)

  implicit none
  integer nrays, nx, ny, nz
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision start_positions(3, nrays), end_positions(3, nrays)
!f2py intent(in) start_positions, end_positions
  real weights(nx,ny,nz)
!f2py intent(in) weights
  real paths(nrays)
!f2py intent(out) paths

  integer n
  real tau


  paths = -999.0
  do n=1,nrays
    tau = 0.0
    !extrapolate ray to domain.
    call util_integrate_ray(start_positions(:,n), end_positions(:,n),&
                            tau, weights, xgrid, ygrid, zgrid,     &
                            nx,ny,nz)
    paths(n) = tau
  enddo
end subroutine util_integrate_rays

subroutine util_integrate_ray(start_position, end_position,       &
  tau, weights, xgrid, ygrid, zgrid, nx,ny,nz)

  implicit none
  integer nrays, nx, ny, nz
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision start_position(3), end_position(3)
!f2py intent(in) start_position, end_position
  real weights(nx,ny,nz), tau, weight1,weight2
!f2py intent(in) weights
!f2py intent(out) tau
  double precision cx, cy, cz, cxinv, cyinv, czinv, eps, r
  double precision sox, soy, soz, so,xn,yn,zn,xe,ye,ze,f(8)
  integer gridptr(3,8), ngrid, maxgrid,i,j
  integer gridptr_end(3,8)
  logical done, a,b,c, final_cell

  tau = 0.0
  eps = 1.0E-5
  ngrid = 0
  maxgrid = 50 * max(nx,ny,nz)

  r = sqrt(sum((end_position - start_position)**2))
  cx = (end_position(1) - start_position(1))/r
  cy = (end_position(2) - start_position(2))/r
  cz = (end_position(3) - start_position(3))/r

  if (abs(cx) .gt. 1.0e-6) then
    cxinv = 1.0D0/cx
  else
    cx = 0.0
    cxinv = 1.0e6
  endif
  if (abs(cy) .gt. 1.0e-6) then
    cyinv = 1.0D0/cy
  else
    cy = 0.0
    cyinv = 1.0e6
  endif
  if (abs(cz) .gt. 1.0e-6) then
    czinv = 1.0D0/cz
  else
    cz = 0.0
    czinv = 1.0e6
  endif

  xe = start_position(1)
  ye = start_position(2)
  ze = start_position(3)

  if (xe .le. xgrid(1)) then
    xe = xgrid(1) + eps
  else if (xe .ge.xgrid(nx)) then
    xe = xgrid(nx) - eps
  endif
  if (ye .le. ygrid(1)) then
    ye = ygrid(1) + eps
  else if (ye .ge.ygrid(ny)) then
    ye = ygrid(ny) - eps
  endif
  if (ze .le. zgrid(1)) then
    ze = zgrid(1) + eps
  else if (ze .ge.zgrid(nz)) then
    ze = zgrid(nz) - eps
  endif

  call util_locate_point(nx, ny, nz, xgrid, ygrid, zgrid,   &
    start_position, gridptr)

  call util_locate_point(nx, ny, nz, xgrid, ygrid, zgrid,   &
    end_position, gridptr_end)

  done = .false.
  do while (.not. done)
    call util_get_interp_kernel(gridptr,xe,ye, ze,f,&
      nx,ny,nz, &
      xgrid, ygrid, zgrid)

    weight1 = 0.0
    do i=1,8
      weight1 = weight1 + f(i)*weights(gridptr(1,i), &
      gridptr(2,i), gridptr(3,i))
    enddo

    if (cx .ge. 1.0e-6) then
      sox = (xgrid(gridptr(1,2)) - xe)*cxinv
    elseif (cx .le. -1.0e-6) then
      sox = (xgrid(gridptr(1,1)) - xe)*cxinv
    else
      sox = 1.0e20
    endif

    if (cy .ge. 1.0e-6) then
      soy = (ygrid(gridptr(2,3)) - ye)*cyinv
    elseif (cy .le. -1.0e-6) then
      soy = (ygrid(gridptr(2,2)) - ye)*cyinv
    else
      soy = 1.0e20
    endif
    if (cz .ge. 1.0e-6) then
      soz = (zgrid(gridptr(3,5)) - ze)*czinv
    elseif (cz .le. -1.0e-6) then
      soz = (zgrid(gridptr(3,4)) - ze)*czinv
    else
      soz = 1.0e20
    endif

    so = min(sox, soy, soz)

    if (so .eq. sox) then
      a = .true.
      b = .false.
      c = .false.
    elseif (so .eq. soy) then
      a = .false.
      b = .true.
      c = .false.
    elseif (so .eq. soz) then
      a = .false.
      b = .false.
      c = .true.
    endif
    xn = xe + so*cx
    yn = ye + so*cy
    zn = ze + so*cz

    final_cell = .true.
    do j=1,8
      do i=1,3
        if (.not. gridptr(i,j) .eq. gridptr_end(i,j)) then
          final_cell = .false.
        endif
      enddo
    enddo

    if (final_cell) then
      xn = end_position(1)
      yn = end_position(2)
      zn = end_position(3)
      done = .true.
    endif

    if (xn > xgrid(nx) + eps) then
      xn = xgrid(nx)
      stop 'went beyond grid +x'
    endif
    if (yn > ygrid(ny) + eps) then
      yn = ygrid(ny)
      stop 'went beyond grid +y'
    endif
    if (zn > zgrid(nz) + eps) then
      zn = zgrid(nz)
      stop 'went beyond grid +z'
    endif
    if (xn < xgrid(1) - eps) then
      xn = xgrid(1)
      stop 'went beyond grid -x'
    endif
    if (yn < ygrid(1)- eps) then
      yn = ygrid(1)
      stop 'went beyond grid -y'
    endif
    if (zn < zgrid(1)- eps) then
      zn = zgrid(1)
      stop 'went beyond grid -z'
    endif

    call util_get_interp_kernel(gridptr,xn,yn, zn,f,&
      nx,ny,nz, &
      xgrid, ygrid, zgrid)

    weight2 = 0.0
    do i=1,8
      weight2 = weight2 + f(i)*weights(gridptr(1,i), &
      gridptr(2,i), gridptr(3,i))
    enddo

    tau = tau + so*0.5*(weight1+weight2)

    if (.not. done) then
      call util_next_cell(a,b,c,cx,cy,cz,gridptr)
    endif

    xe = xn
    ye = yn
    ze = zn
    ngrid = ngrid + 1
    if (ngrid .ge. maxgrid) then
      done = .true.
    endif

  enddo

end subroutine util_integrate_ray

subroutine util_next_cell(a,b,c,cx,cy,cz,gridptr)

  implicit none
  logical a,b,c
!f2py intent(in) a,b,c
  double precision cx, cy, cz
!f2py intent(in) cx, cy, cz
  integer gridptr(3,8)
!f2py intent(in, out) gridptr

  integer i

  if (a) then
    if (cx .ge. 1e-6) then
      do i=1,8
        gridptr(1,i) = gridptr(1,i) + 1
      enddo
    elseif (cx .le. -1e-6) then
      do i=1,8
        gridptr(1,i) = gridptr(1,i) - 1
      enddo
    endif
  endif
  if (b) then
    if (cy .ge. 1e-6) then
      do i=1,8
        gridptr(2,i) = gridptr(2,i) + 1
      enddo
    elseif (cy .le. -1e-6) then
      do i=1,8
        gridptr(2,i) = gridptr(2,i) - 1
      enddo
    endif
  endif
  if (c) then
    if (cz .ge. 1e-6) then
      do i=1,8
        gridptr(3,i) = gridptr(3,i) + 1
      enddo
    elseif (cz .le. -1e-6) then
      do i=1,8
        gridptr(3,i) = gridptr(3,i) - 1
      enddo
    endif
  endif
  return

end subroutine util_next_cell

! subroutine extrapolate_to_domain_edge(nx,ny,nz,xgrid,ygrid,zgrid, &
!   cx_old,cy_old,cz_old,cxinv,cyinv,czinv, position, end)
!   implicit none
!
!   integer nx,ny,nz
!   double precision xgrid(nx), ygrid(ny), zgrid(nz)
! !f2py intent(in) xgrid, ygrid, zgrid
!   double precision position(3),cx_old,cy_old,cz_old
!   double precision cxinv, cyinv, czinv
! !f2py intent(in, out) position
! !f2py intent(in) cx_old, cy_old, cz_old
! !f2py intent(in) cxinv, cyinv, czinv
!   logical end
! !f2py intent(in) end
!
!   double precision sox1,sox2,soy1,soy2,soz1,soz2,so
!   double precision cx,cy,cz,x0,y0,z0
!
! !if outside domain, extrapolate ray to domain edge.
!   if ((xgrid(1) .ge. position(1)) .or. &
!     (xgrid(nx) .le. position(1)) .or. &
!     (ygrid(1) .ge. position(2)) .or. &
!       (ygrid(ny) .le. position(2)) .or. &
!     (zgrid(1) .ge. position(3)) .or. &
!       (zgrid(nz) .le. position(3))) then
!
!     if (end) then
!       cx = -cx_old
!       cy = -cy_old
!       cz = -cz_old
!     else
!       cx = cx_old
!       cy = cy_old
!       cz = cz_old
!     endif
!
! !distances to all potential matches.
!     if (cx .ge. 1.0e-6) then
!       sox1 = (xgrid(1) - position(1))*cxinv
!       sox2 = (xgrid(nx) - position(1))*cxinv
!     elseif (cx .le. -1.0e-6) then
!       sox1 = (xgrid(1) - position(1))*cxinv
!       sox2 = (xgrid(nx) - position(1))*cxinv
!     else
!       sox1 = 1.0e20
!       sox2 = 1.0e20
!     endif
!     if (cy .ge. 1.0e-6) then
!       soy1 = (ygrid(1) - position(2))*cyinv
!       soy2= (ygrid(ny) - position(2))*cyinv
!     elseif (cy .le. -1.0e-6) then
!       soy1 = (ygrid(1) - position(2))*cyinv
!       soy2 = (ygrid(ny) - position(2))*cyinv
!     else
!       soy1 = 1.0e20
!       soy2 = 1.0e20
!     endif
!     if (cz .ge. 1.0e-6) then
!       soz1 = (zgrid(1) - position(3))*czinv
!       soz2 = (zgrid(nz) - position(3))*czinv
!     elseif (cz .le. -1.0e-6) then
!       soz1 = (zgrid(1) - position(3))*czinv
!       soz2 = (zgrid(nz) - position(3))*czinv
!     else
!       soz1 = 1.0e20
!       soz2 = 1.0e20
!     endif
!     print *, sox1,sox2,soy1,soy2,soz1,soz2
! !test each potential match for validity.
! !take the valid match with the smallest distance.
!
!     so = 1.0e20
!     if (min(sox1, so) .eq. sox1) then
!       z0 = position(3) + sox1*cz
!       y0 = position(2) + sox1*cy
!       if ((z0 .ge. zgrid(1)) .and. (z0 .le. zgrid(nz)) .and. &
!         (y0 .ge. ygrid(1)) .and. (y0 .le. ygrid(ny))) then
!         so = sox1
!       endif
!     endif
!       print *,'sox1', sox1,so
!     if (min(sox2, so) .eq. sox2) then
!       z0 = position(3) + sox2*cz
!       y0 = position(2) + sox2*cy
!       if ((z0 .ge. zgrid(1)) .and. (z0 .le. zgrid(nz)) .and. &
!         (y0 .ge. ygrid(1)) .and. (y0 .le. ygrid(ny))) then
!         so = sox2
!       endif
!     endif
!       print *,'sox2', sox2,so
!     if (min(soy1, so) .eq. soy1) then
!       z0 = position(3) + soy1*cz
!       x0 = position(1) + soy1*cx
!       if ((z0 .ge. zgrid(1)) .and. (z0 .le. zgrid(nz)) .and. &
!         (x0 .ge. xgrid(1)) .and. (x0 .le. xgrid(nx))) then
!         so = soy1
!       endif
!     endif
!       print *,'soy1', soy1,so
!     if (min(soy2, so) .eq. soy2) then
!       z0 = position(3) + soy2*cz
!       x0 = position(1) + soy2*cx
!       if ((z0 .ge. zgrid(1)) .and. (z0 .le. zgrid(nz)) .and. &
!         (x0 .ge. xgrid(1)) .and. (x0 .le. xgrid(nx))) then
!         so = soy2
!       endif
!     endif
!       print *,'soy2', soy2,so
!     if (min(soz1, so) .eq. soz1) then
!       y0 = position(2) + soz1*cy
!       x0 = position(1) + soz1*cx
!       if ((y0 .ge. ygrid(1)) .and. (y0 .le. ygrid(ny)) .and. &
!         (x0 .ge. xgrid(1)) .and. (x0 .le. xgrid(nx))) then
!         so = soz1
!       endif
!     endif
!       print *,'soz1', soz1,so
!     if (min(soz2, so) .eq. soz2) then
!       y0 = position(2) + soz2*cy
!       x0 = position(1) + soz2*cx
!       if ((y0 .ge. ygrid(1)) .and. (y0 .le. ygrid(ny)) .and. &
!         (x0 .ge. xgrid(1)) .and. (x0 .le. xgrid(nx))) then
!         so = soy2
!       endif
!     endif
!       print *,'soz2', soz2,so
!
!     position(1) = position(1) + so*cx
!     position(2) = position(2) + so*cy
!     position(3) = position(3) + so*cz
!
!   endif
!   return
! end subroutine extrapolate_to_domain_edge
!

subroutine util_locate_point(nx, ny, nz, xgrid, ygrid, zgrid, &
  location, gridptr)
  implicit none

  integer nx, ny, nz, gridptr(3,8)
!f2py intent(out) gridptr
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision location(3)
!f2py intent(in) location

  integer i,j
  gridptr = 1

  call nearest_binary(nx, xgrid, location(1),i,j)
  gridptr(1,1) = i
  gridptr(1,2) = j
  gridptr(1,3) = i
  gridptr(1,4) = j
  gridptr(1,5) = i
  gridptr(1,6) = j
  gridptr(1,7) = i
  gridptr(1,8) = j

  call nearest_binary(ny, ygrid, location(2),i,j)
  gridptr(2,1) = i
  gridptr(2,2) = i
  gridptr(2,3) = j
  gridptr(2,4) = j
  gridptr(2,5) = i
  gridptr(2,6) = i
  gridptr(2,7) = j
  gridptr(2,8) = j

  call nearest_binary(nz, zgrid, location(3),i,j)
  gridptr(3,1) = i
  gridptr(3,2) = i
  gridptr(3,3) = i
  gridptr(3,4) = i
  gridptr(3,5) = j
  gridptr(3,6) = j
  gridptr(3,7) = j
  gridptr(3,8) = j

  ! do i=1,nx
  !   if ((location(1) .ge. xgrid(i)) .and.  &
  !         (location(1) .le. xgrid(i+1))) then
  !     gridptr(1,1) = i
  !     gridptr(1,2) = i + 1
  !     gridptr(1,3) = i + 1
  !     gridptr(1,4) = i
  !     gridptr(1,5) = i
  !     gridptr(1,6) = i + 1
  !     gridptr(1,7) = i + 1
  !     gridptr(1,8) = i
  !   endif
  ! enddo
  !
  ! do i=1,ny
  !   if ((location(2) .ge. ygrid(i)) .and.  &
  !         (location(2) .le. ygrid(i+1))) then
  !     gridptr(2,1) = i
  !     gridptr(2,2) = i
  !     gridptr(2,3) = i + 1
  !     gridptr(2,4) = i + 1
  !     gridptr(2,5) = i
  !     gridptr(2,6) = i + 1
  !     gridptr(2,7) = i + 1
  !     gridptr(2,8) = i
  !   endif
  ! enddo
  !
  ! do i=1,nz
  !   if ((location(3) .ge. zgrid(i)) .and.  &
  !         (location(3) .le. zgrid(i+1))) then
  !     gridptr(3,1) = i
  !     gridptr(3,2) = i
  !     gridptr(3,3) = i
  !     gridptr(3,4) = i
  !     gridptr(3,5) = i + 1
  !     gridptr(3,6) = i + 1
  !     gridptr(3,7) = i + 1
  !     gridptr(3,8) = i + 1
  !   endif
  ! enddo
  return
end subroutine util_locate_point

subroutine util_get_interp_kernel2(i,j,k,xn,yn, zn,f,nx,ny,nz, &
  xgrid, ygrid, zgrid)
  implicit none

  integer nx,ny,nz,i,j,k!index_low(3), index_high(3)
!f2py intent(in) i,j,k
  double precision xn, yn, zn, f(8)
!f2py intent(in) xn, yn, zn
!f2py intent(out) f
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision u,v,w, delx,dely, delz
  double precision invdelx, invdely, invdelz

  if ((xn > xgrid(i+1)) .or. (xn < xgrid(i)) .or. &
      (yn > ygrid(j+1)) .or. (yn < ygrid(j)) .or. &
      (zn > zgrid(k+1)) .or. (zn < zgrid(k))) then
    f = 0.0D0
  else

    delx = xgrid(i + 1) - xgrid(i)!xgrid(index_high(1)) - xgrid(index_low(1))
    dely = ygrid(j + 1) - ygrid(j)!ygrid(index_high(2)) - xgrid(index_low(2))
    delz = zgrid(k + 1) - zgrid(k)!zgrid(index_high(3)) - xgrid(index_low(3))

    if (delx .le. 0.0) then
      invdelx = 1.0
    else
      invdelx = 1.0D0/delx
    endif
    if (dely .le. 0.0) then
      invdely = 1.0
    else
      invdely = 1.0D0/dely
    endif
    if (delz .le. 0.0) then
      invdelz = 1.0
    else
      invdelz = 1.0D0/delz
    endif

    u = (xn - xgrid(i))*invdelx
    v = (yn - ygrid(j))*invdely
    w = (zn - zgrid(k))*invdelz

    f(1) = (1-w)*(1-v)*(1-u)
    f(2) = (1-w)*(1-v)*u
    f(3) = (1-w)*v*(1-u)
    f(4) = (1-w)*v*u
    f(5) = w*(1-v)*(1-u)
    f(6) = w*(1-v)*u
    f(7) = w*v*(1-u)
    f(8) = w*v*u
  endif
  return
end subroutine util_get_interp_kernel2




subroutine util_get_interp_kernel(gridptr,xn,yn, zn,f,nx,ny,nz, &
  xgrid, ygrid, zgrid)
  implicit none

  integer nx,ny,nz, gridptr(3,8)
!f2py intent(in) gridptr
  double precision xn, yn, zn, f(8)
!f2py intent(in) xn, yn, zn
!f2py intent(out) f
  double precision xgrid(nx), ygrid(ny), zgrid(nz)
!f2py intent(in) xgrid, ygrid, zgrid
  double precision u,v,w, delx,dely, delz
  double precision invdelx, invdely, invdelz

  delx = xgrid(gridptr(1,2)) - xgrid(gridptr(1,1))
  dely = ygrid(gridptr(2,3)) - xgrid(gridptr(2,2))
  delz = zgrid(gridptr(3,5)) - xgrid(gridptr(3,4))

  if (delx .le. 0.0) then
    invdelx = 1.0
  else
    invdelx = 1.0D0/delx
  endif
  if (dely .le. 0.0) then
    invdely = 1.0
  else
    invdely = 1.0D0/dely
  endif
  if (delz .le. 0.0) then
    invdelz = 1.0
  else
    invdelz = 1.0D0/delz
  endif

  u = (xn - xgrid(gridptr(1,1)))*invdelx
  v = (yn - ygrid(gridptr(2,1)))*invdely
  w = (zn - zgrid(gridptr(3,1)))*invdelz

  f(1) = (1-W)*(1-V)*(1-U)
  f(2) = (1-W)*(1-V)*U
  f(3) = (1-W)*V*(1-U)
  f(4) = (1-W)*V*U
  f(5) = W*(1-V)*(1-U)
  f(6) = W*(1-V)*U
  f(7) = W*V*(1-U)
  f(8) = W*V*U

  return
end subroutine util_get_interp_kernel


subroutine test_source(NSTOKES, NPTS, ML, MM, NLM, NLEG, NSTLEG, &
  DELTAM, NUMPHASE, NPART, MAXNMICRO,IPHASE, PHASEINTERPWT, &
  INTERPMETHOD, PHASEMAX,YLMSUN, SOLARMU,ALBEDO, LEGEN, EXTINCT, &
  TOTAL_EXT, PLANCK, DIRFLUX, RADIANCE, SOURCE, SRCTYPE, IERR, &
  ERRMSG, NEWMETHOD, TIME_OUT)
  IMPLICIT NONE
  INTEGER NSTOKES, NPTS, ML, MM, NLM, NLEG, NSTLEG
!f2py intent(in) :: NSTOKES, NPTS, ML, MM, NLM, NLEG, NSTLEG
  LOGICAL DELTAM
!f2py intent(in) :: DELTAM
  INTEGER NUMPHASE, NPART, MAXNMICRO
!f2py intent(in) :: NUMPHASE, NPART, MAXNMICRO
  INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART)
  REAL    PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART)
  CHARACTER INTERPMETHOD*2
  REAL    PHASEMAX
!f2py intent(in) :: IPHASE, PHASEINTERPWT, INTERPMETHOD, PHASEMAX
  REAL    YLMSUN(NSTLEG,NLM), SOLARMU
!f2py intent(in) :: YLMSUN, SOLARMU
  REAL    ALBEDO(NPTS,NPART), LEGEN(NSTLEG,0:NLEG,*)
  REAL    EXTINCT(NPTS,NPART), TOTAL_EXT(*)
!f2py intent(in) :: ALBEDO, LEGEN, EXTINCT, TOTAL_EXT
  REAL    PLANCK(NPTS,NPART), DIRFLUX(*)
!f2py intent(in) :: PLANCK, DIRFLUX
  REAL    RADIANCE(NSTOKES,NLM)
  REAL    SOURCE(NSTOKES,NLM)
!f2py intent(in) :: RADIANCE
!f2py intent(out) :: SOURCE
  CHARACTER SRCTYPE*1
!f2py intent(in) :: SRCTYPE
  INTEGER IERR
  CHARACTER ERRMSG*600
!f2py intent(out) :: IERR, ERRMSG
  LOGICAL NEWMETHOD
!f2py intent(in) :: NEWMETHOD
  DOUBLE PRECISION TIME_OUT
!f2py intent(out) :: TIME_OUT

  INTEGER IS, ISO, IR, I, IPH, J, JS, K, L, LS, M, MS, ME, NS, NR
  INTEGER IPA, Q
  DOUBLE PRECISION TIME1, TIME2
  INTEGER, ALLOCATABLE :: LOFJ(:)
  REAL, ALLOCATABLE :: SOURCET(:,:),  SOURCET1(:,:)
  REAL    SRCMIN, C, SECMU0, D, EXT, W, F, ALB, SCAT, TOTAL_PLANCK
  REAL, ALLOCATABLE :: LEGENT(:,:,:), LEGENT1(:,:,:)
  ALLOCATE (LEGENT(NSTLEG,0:NLEG,1),LEGENT1(NSTLEG,0:NLEG,1))
  ALLOCATE (SOURCET(NSTOKES, NLM), LOFJ(NLM), SOURCET1(NSTOKES,NLM))

  CALL CPU_TIME(TIME1)

  SECMU0 = 1.0D0/ABS(SOLARMU)
  I = 1

  J = 0
  DO L = 0, ML
    ME = MIN(L,MM)
    MS = -ME
    DO M = MS, ME
      J = J + 1
      LOFJ(J) = L
    ENDDO
  ENDDO

  IF (NEWMETHOD) THEN
    EXT = TOTAL_EXT(I)
    ALB = 0.0
    LEGENT = 0.0
    TOTAL_PLANCK = 0.0

    DO IPA = 1, NPART
     SCAT = EXTINCT(I,IPA)*ALBEDO(I,IPA)
     ALB = ALB + SCAT
     TOTAL_PLANCK = TOTAL_PLANCK + &
        EXTINCT(I,IPA)*PLANCK(I,IPA)
     IF (INTERPMETHOD(2:2) .EQ. 'O') THEN
       LEGENT(:,:,1) = LEGENT(:,:,1) + &
          SCAT*LEGEN(:,:,IPHASE(1,I,IPA))
     ELSEIF (INTERPMETHOD(2:2) .EQ. 'N') THEN
       IF (PHASEINTERPWT(1,I,IPA) .GE. PHASEMAX) THEN
         LEGENT(:,:,1) = LEGENT(:,:,1) + &
            SCAT*LEGEN(:,:,IPHASE(1,I,IPA))
       ELSE
         LEGENT1(:,:,:) = 0.0
         DO Q=1,8*MAXNMICRO
           IF (PHASEINTERPWT(Q,I,IPA) .LE. 1e-5) CYCLE
           IPH = IPHASE(Q,I,IPA)
           LEGENT1(:,:,1) = LEGENT1(:,:,1) + &
            LEGEN(:,:,IPH)* &
            PHASEINTERPWT(Q,I,IPA)
         ENDDO
         LEGENT(:,:,:) = LEGENT(:,:,:) + &
            SCAT*LEGENT1(:,:,:)
       ENDIF
     ENDIF
    ENDDO
    IF (ALB .NE. 0.0) THEN
     LEGENT(:,:,:) = LEGENT(:,:,:)/ALB
    ENDIF
    IF (INTERPMETHOD(2:2) .EQ. 'N') THEN
      IF (DELTAM) THEN
        F = LEGENT(1,ML+1,1)
        LEGENT(:,0:ML,1) = LEGENT(:,0:ML,1)/(1-F)
      ENDIF
    ENDIF

    IF (EXT .GT. 1e-10) THEN
      ALB = ALB/EXT
      TOTAL_PLANCk = TOTAL_PLANCK/EXT
    ELSE
      ALB = 0.0
      TOTAL_PLANCK = 0.0
    ENDIF

    IF (NSTOKES .EQ. 1) THEN
      CALL CALC_SOURCE_PNT_UNPOL (NLM, NLEG, LOFJ, &
          SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN, &
          TOTAL_PLANCK, ALB, 1, LEGENT,&
          NLM, RADIANCE(1,1),   SOURCET)
    ELSE
      CALL CALC_SOURCE_PNT (NSTOKES, NLM, NSTLEG, NLEG,&
          LOFJ, SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
          TOTAL_PLANCK, ALB, 1, LEGENT,&
          NLM, RADIANCE(1,1), SOURCET)
    ENDIF

    ! PRINT *, SOURCET

  ELSE

    EXT = TOTAL_EXT(I)
    SOURCET = 0.0
    DO IPA = 1, NPART
      IF (EXT.EQ.0.0) THEN
        W = 1.0
      ELSE
        W = EXTINCT(I,IPA)/EXT
      ENDIF
      IF (INTERPMETHOD(2:2) .EQ. 'O') THEN
        IPH = IPHASE(1,I,IPA)
        IF (NSTOKES .EQ. 1) THEN
          CALL CALC_SOURCE_PNT_UNPOL (NLM, NLEG, LOFJ, &
            SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN, &
            PLANCK(I,IPA), ALBEDO(I,IPA), IPH, LEGEN, &
            NLM, RADIANCE(1,1),   SOURCET1)
        ELSE
          CALL CALC_SOURCE_PNT (NSTOKES, NLM, NSTLEG, NLEG, &
            LOFJ, SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN, &
            PLANCK(I,IPA), ALBEDO(I,IPA), IPH, LEGEN, &
            NLM, RADIANCE(1,1), SOURCET1)
        ENDIF
      ELSEIF (INTERPMETHOD(2:2) .EQ. 'N') THEN
        IF (PHASEINTERPWT(1,I,IPA) .GE. PHASEMAX) THEN
          LEGENT(:,:,1)  = LEGEN(:,:,IPHASE(1,I,IPA))
        ELSE
          LEGENT = 0.0
          DO Q=1,8*MAXNMICRO
            IF (PHASEINTERPWT(Q,I,IPA) .LE. 1e-5) CYCLE
            IPH = IPHASE(Q,I,IPA)
            LEGENT(:,:,1) = LEGENT(:,:,1) + &
            LEGEN(:,:,IPH)* &
            PHASEINTERPWT(Q,I,IPA)
          ENDDO
        ENDIF
        IF (DELTAM) THEN
          F = LEGENT(1,ML+1,1)
          LEGENT(:,0:ML,1) = LEGENT(:,0:ML,1)/(1-F)
        ENDIF
        IF (NSTOKES .EQ. 1) THEN
          CALL CALC_SOURCE_PNT_UNPOL (NLM, NLEG, LOFJ, &
            SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN, &
            PLANCK(I,IPA), ALBEDO(I,IPA), 1, LEGENT, &
            NLM, RADIANCE(1,1),   SOURCET1)
        ELSE
          CALL CALC_SOURCE_PNT (NSTOKES, NLM, NSTLEG, NLEG, &
            LOFJ, SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN, &
            PLANCK(I,IPA), ALBEDO(I,IPA), 1, LEGENT, &
            NLM, RADIANCE(1,1), SOURCET1)
        ENDIF
      ENDIF
      SOURCET = SOURCET + W * SOURCET1
    ENDDO

  ENDIF
  SOURCE(:,:) = SOURCET(:,:)
  DEALLOCATE(SOURCET, SOURCET1, LEGENT, LEGENT1, LOFJ)
  CALL CPU_TIME(TIME2)
  TIME_OUT = DBLE(TIME2) - DBLE(TIME1)
END

! A light integrator for emission problems.
! Useful for atmospheric correction problems.
! Assume that the (1) index is the TOA.
subroutine integrate_thermal_source(extinction, planck,z,&
                            tautol, cumulative_radiance, &
                            nz, nprof, surface_emission, &
                            cumulative_transmit)
  implicit none
  integer :: nz, nprof
!f2py intent(in) :: nz, nprof
  real*8 :: z(nz, nprof), extinction(nz, nprof), planck(nz, nprof)
!f2py intent(in) :: z, extinction, planck
  real*8 :: cumulative_radiance(nz, nprof), surface_emission(nprof)
!f2py intent(in) :: surface_emission
!f2py intent(out) :: cumulative_radiance
  real*8 :: cumulative_transmit(nz, nprof)
!f2py intent(out) :: cumulative_transmit
  real*8 :: tautol
!f2py intent(in) :: tautol

  integer i, iz, iprof, ntau

  real*8 :: inv_delz(nz - 1), diff, so, taugrid, dels, ext1, srcext1
  real*8 :: s, u, ext0,srcext0, tau, abscell, transcell, src, ext, transmit

  cumulative_radiance = 0.0d0

  do iprof = 1, nprof
    transmit = 1.0d0
    cumulative_transmit(1,iprof) = 1.0d0

    do iz=1,nz-1
      diff = abs(z(iz,iprof) - z(iz + 1,iprof))
      if (diff .gt. 1e-16) then
        inv_delz(iz) = 1.0d0/diff
      else
        inv_delz(iz) = 1e16
      endif
    enddo

    do iz = 1,nz - 1 ! cell indices
      so = abs(z(iz,iprof) - z(iz + 1,iprof))
      taugrid = so*0.5d0*(extinction(iz + 1,iprof)+extinction(iz,iprof))
      ntau = max(1,1+int(taugrid/tautol))
      dels = so/ntau

      ext1 = extinction(iz,iprof)
      srcext1 = ext1*planck(iz,iprof)

!      print *, iprof, iz, ext1, planck(iz,iprof)

      do i=1,ntau
        s = i*dels
        u = max(min(s*inv_delz(iz),1.0d0),0.0d0)
        ext0 = u*extinction(iz + 1,iprof) + (1.0d0-u)*extinction(iz,iprof)
        srcext0 = u*extinction(iz + 1,iprof)*planck(iz + 1,iprof) + &
                  (1.0d0-u)*extinction(iz, iprof)*planck(iz,iprof)
        ext = 0.5d0*(ext0 + ext1)
        if (ext .gt. 0.0d0) then
          tau = ext*dels
          abscell = tau*(1.0d0 - 0.5d0*tau*(1.0d0-0.333333333333d0*tau))
          transcell = 1.0d0 - abscell
          src = ( 0.5d0*(srcext0 + srcext1) + &
            0.08333333333d0*(ext0*srcext1-ext1*srcext0)*&
            dels*(1.0d0 - 0.05d0 *(ext1-ext0)*dels))/ ext
        else
          abscell = 0.0d0
          transcell = 1.0d0
          src = 0.0d0
        endif
!        print *, iprof, iz, i, ext, tau, dels, src, abscell, transcell, transmit
        cumulative_radiance(iz + 1,iprof) = cumulative_radiance(iz + 1,iprof) +&
                                           transmit*src*abscell
        transmit = transmit*transcell
        cumulative_transmit(iz + 1,iprof) = transmit
        ext1 = ext0
        srcext1 = srcext0
      enddo
    enddo
    cumulative_radiance(nz,iprof) = cumulative_radiance(nz,iprof) + &
                transmit*surface_emission(iprof)
  enddo
end subroutine integrate_thermal_source